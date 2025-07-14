#include "problem/problem_proxpde.hpp"

// std
#include <string_view>

// proxpde
#include <proxpde/assembly.hpp>
#include <proxpde/builder.hpp>

// local
#include "enums.hpp"

namespace cocoa
{

std::unique_ptr<Problem> ProblemProXPDE::build(EQN_TYPE const type)
{
  switch (type)
  {
  case EQN_TYPE::HEAT:
  case EQN_TYPE::HEAT_COUPLED:
  case EQN_TYPE::HEAT_BUOYANT:
    return std::unique_ptr<Problem>(new ProblemProXPDEHeat);
  case EQN_TYPE::NS:
  case EQN_TYPE::NS_BOUSSINESQ:
    return std::unique_ptr<Problem>(new ProblemProXPDENS);
  default:
    std::abort();
  }
}

ProblemProXPDE::ProblemProXPDE(): Problem{PROBLEM_TYPE::PROXPDE} {}

std::unique_ptr<MeshCoupling> ProblemProXPDE::initMeshCoupling(
    COUPLING_TYPE type, COUPLING_SCOPE scope, Marker marker, std::string_view bdName)
{
  if (scope == COUPLING_SCOPE::VOLUME)
  {
    using Elem_T = Mesh_T::Elem_T;

    // coords format: x_0, y_0, z_0, x_1, ...
    std::vector<double> coords(mesh_.pointList.size() * 3);
    for (auto const & p: mesh_.pointList)
    {
      coords[3 * p.id] = p.coord[0];
      coords[3 * p.id + 1] = p.coord[1];
      coords[3 * p.id + 2] = p.coord[2];
    }

    // conn format: elem0_numpts, id_0, id_1, ..., elem1_numpts, ...
    std::vector<uint> conn(mesh_.elementList.size() * (Elem_T::numPts + 1));
    for (auto const & elem: mesh_.elementList)
    {
      conn[(Elem_T::numPts + 1) * elem.id] = MEDCellTypeToIKCell(MED_CELL_TYPE::QUAD4);
      for (uint p = 0; p < Elem_T::numPts; p++)
      {
        conn[(Elem_T::numPts + 1) * elem.id + p + 1] = elem.pts[p]->id;
      }
    }

    // offsets format: sum_0^k elemk_numpts + 1,
    std::vector<uint> offsets(mesh_.elementList.size() + 1);
    offsets[0] = 0;
    for (uint e = 0; e < mesh_.elementList.size(); e++)
    {
      offsets[e + 1] = offsets[e] + (Elem_T::numPts + 1);
    }

    auto meshCoupling = MeshCoupling::build(type);
    meshCoupling->init(
        "mesh_proxpde",
        COUPLING_SCOPE::VOLUME,
        markerNotSet,
        "",
        Elem_T::dim,
        3u,
        coords,
        conn,
        offsets);
    return meshCoupling;
  }
  else if (scope == COUPLING_SCOPE::BOUNDARY)
  {
    fmt::println(stderr, "boundary coupling not yet implemented!");
    std::abort();
  }
  else
  {
    std::abort();
  }
}

void ProblemProXPDE::advance()
{
  if (time + dt_ < finalTime_ - 1.e-6)
  {
    time += dt_;
  }
  else
  {
    dt_ = finalTime_ - time;
    time = finalTime_;
  }
  it++;
}

bool ProblemProXPDE::run() const { return time < finalTime_; }

uint ProblemProXPDE::solve()
{
  fmt::println("{} - time = {:.6e}, dt = {:.6e}", name_, time, dt_);

  // update old values
  uOld_ = u_;

  // build matrix and rhs
  proxpde::Builder<> builder{this->size()};
  assemblies_.at(equationType_)->evaluate(this, builder);

  // solve
  proxpde::LUSolver solver;
  solver.compute(builder.A);
  u_ = solver.solve(builder.b);

  // TODO: extract number of iterations
  return 0u;
}

Marker ProblemProXPDE::findRegion(std::string_view name)
{
  if (name == "left")
    return 0u;
  else if (name == "right")
    return 1u;
  else if (name == "bottom")
    return 2u;
  else if (name == "top")
    return 3u;
  else
  {
    fmt::println("region {} not recognized", name);
    std::abort();
  }
  return markerNotSet;
}

// =====================================================================

struct AssemblyHeat: public ProblemProXPDE::Assembly
{
  auto evaluate(ProblemProXPDE * pParent, proxpde::Builder<> & b) -> void override
  {
    auto p = dynamic_cast<ProblemProXPDEHeat *>(pParent);
    // update coupling
    // setData("vel", p->vel_.data, *p->vel_.feSpace);
    auto const alpha = p->params_.at("alpha");
    auto const q = p->fieldsP0_.at("q");

    // assembly lhs
    b.buildLhs(
        std::tuple{
            proxpde::AssemblyMass{1. / p->dt_, p->feSpace_},
            proxpde::AssemblyStiffness{alpha, p->feSpace_},
            proxpde::AssemblyAdvection{p->vel_, p->feSpace_},
        },
        p->bcs_);
    // fmt::println(stderr, "triplets");
    // for (auto const & t: b._triplets)
    // {
    //   fmt::println(stderr, "{}, {} -> {}", t.row(), t.col(), t.value());
    // }
    b.closeMatrix();

    // assembly rhs
    b.buildRhs(
        std::tuple{
            proxpde::AssemblyProjection{1. / p->dt_, p->uOld_, p->feSpace_},
            proxpde::AssemblyProjection{1., q.data, p->feSpaceP0_, p->feSpace_},
        },
        p->bcs_);
  }
};

struct AssemblyHeatCoupled: public ProblemProXPDE::Assembly
{
  auto evaluate(ProblemProXPDE * pParent, proxpde::Builder<> & b) -> void override
  {
    auto p = dynamic_cast<ProblemProXPDEHeat *>(pParent);
    auto const alpha = p->params_.at("alpha");
    auto const q = p->fieldsP0_.at("q");

    // // update coupling
    // p->getDataMED("vel", p->vel_.data, *p->vel_.feSpace);

    // lhs terms
    auto timeDer = proxpde::AssemblyMass{1. / p->dt_, p->feSpace_};
    auto diffusion = proxpde::AssemblyStiffness{alpha, p->feSpace_};
    auto advection = proxpde::AssemblyAdvection{p->vel_, p->feSpace_};

    proxpde::Vec mask(p->feSpace_.dof.size);
    // p->getDataMED("mask", mask, p->feSpace_);
    auto const kAmpli = 20.;
    proxpde::FEVar kAmpliVec("kAmpli", p->feSpace_);
    kAmpliVec.data = kAmpli * mask;
    auto feedbackLhs = proxpde::AssemblyMassFE{1.0, kAmpliVec, p->feSpace_};

    // rhs terms
    auto timeDerOld = proxpde::AssemblyProjection{1. / p->dt_, p->uOld_, p->feSpace_};
    auto source = proxpde::AssemblyProjection{1., q.data, p->feSpaceP0_, p->feSpace_};

    proxpde::Vec Tcfd(p->feSpace_.dof.size);
    // p->getDataMED("T", Tcfd, p->feSpace_);
    Tcfd.array() *= mask.array();
    auto feedbackRhs = proxpde::AssemblyProjection{kAmpli, Tcfd, p->feSpace_};

    // assembly
    b.buildLhs(std::tuple{timeDer, diffusion, advection, feedbackLhs}, p->bcs_);
    b.closeMatrix();
    b.buildRhs(std::tuple{timeDerOld, source, feedbackRhs}, p->bcs_);
  }
};

struct AssemblyHeatBuoyant: public ProblemProXPDE::Assembly
{
  auto evaluate(ProblemProXPDE * pParent, proxpde::Builder<> & b) -> void override
  {
    auto p = dynamic_cast<ProblemProXPDEHeat *>(pParent);
    auto const alpha = p->params_.at("alpha");

    // update coupling
    // p->getDataMED("vel", p->vel_.data, *p->vel_.feSpace);
    // p->vel_ << [](proxpde::Vec3 const & pt) {
    //   return proxpde::Vec2{pt[1] - 0.5, 0.5 - pt[0]};
    // };
    // p->vel_ << [](proxpde::Vec3 const & pt)
    // {
    //   return proxpde::Vec2{
    //       -2.0 * std::sin(M_PI * pt[0]) * std::sin(M_PI * pt[0]) *
    //           std::sin(M_PI * pt[1]) * std::cos(M_PI * pt[1]),
    //       2.0 * std::sin(M_PI * pt[0]) * std::cos(M_PI * pt[0]) *
    //           std::sin(M_PI * pt[1]) * std::sin(M_PI * pt[1]),
    //   };
    // };

    // lhs terms
    auto timeDer = proxpde::AssemblyMass{1. / p->dt_, p->feSpace_};
    auto diffusion = proxpde::AssemblyStiffness{alpha, p->feSpace_};
    auto advection = proxpde::AssemblyAdvection{p->vel_, p->feSpace_};

    // rhs terms
    auto timeDerOld = proxpde::AssemblyProjection{1. / p->dt_, p->uOld_, p->feSpace_};
    // auto source = proxpde::AssemblyProjection{1., q.data, p->feSpace_};

    // assembly lhs
    b.buildLhs(std::tuple{timeDer, diffusion, advection}, p->bcs_);
    b.closeMatrix();

    // assembly rhs
    b.buildRhs(
        std::tuple{
            timeDerOld,
            // source,
        },
        p->bcs_);
  }
};

ProblemProXPDEHeat::ProblemProXPDEHeat()
{
  [[maybe_unused]] auto const & [itHeat, successHeat] =
      assemblies_.emplace(EQN_TYPE::HEAT, new AssemblyHeat);
  assert(successHeat);

  [[maybe_unused]] auto const & [itHC, successHC] =
      assemblies_.emplace(EQN_TYPE::HEAT_COUPLED, new AssemblyHeatCoupled);
  assert(successHC);

  [[maybe_unused]] auto const & [itHB, successHB] =
      assemblies_.emplace(EQN_TYPE::HEAT_BUOYANT, new AssemblyHeatBuoyant);
  assert(successHB);
}

void ProblemProXPDEHeat::setup(Problem::ConfigList_T const & configs)
{
  // get configuration from file
  // TODO: string conversion is required for python bindings, manage other variant types
  // as errors
  auto const & configFileVariant = configs.at("config_file");
  auto const configFile =
      std::holds_alternative<std::filesystem::path>(configFileVariant)
          ? std::get<std::filesystem::path>(configFileVariant).string()
          : std::get<std::string>(configFileVariant);
  proxpde::ParameterDict config = YAML::LoadFile(configFile);

  name_ = config["name"].as<std::string>();

  // init mesh from configuration
  proxpde::readMesh(mesh_, config["mesh"]);

  // init equation
  equationType_ = str2eqn(config["equation_type"].as<std::string>());

  // init time related members
  time = 0.0;
  dt_ = config["dt"].as<double>();
  finalTime_ = config["final_time"].as<double>();

  // init fe related members
  feSpace_.init(mesh_);
  feSpaceP0_.init(mesh_);
  T_.init("T", feSpace_);
  double const initValue = config["init_value"].as<double>();
  T_ << initValue;
  u_ = T_.data;

  // init physical constants and operating conditions
  if (config["alpha"])
    params_["alpha"] = config["alpha"].as<double>();
  if (config["q"])
  {
    fieldsP0_.emplace("q", proxpde::FEVar<FESpaceP0_T>{"q", feSpaceP0_});
    fieldsP0_.at("q") << config["q"].as<double>();
  }
  feSpaceVel_.init(mesh_);
  vel_.init("vel", feSpaceVel_);
  vel_ << config["vel"].as<proxpde::Vec2>();

  // init additional fields
  if (config["additional_fields_p1"])
  {
    for (auto const & name:
         config["additional_fields_p1"].as<std::vector<std::string>>())
    {
      auto const & [it, success] =
          fieldsP1_.emplace(name, proxpde::FEVar<FESpace_T>{name, feSpace_});
      assert(success);
      it->second.data = proxpde::Vec::Zero(feSpace_.dof.size);
    }
  }
  if (config["additional_fields_p0"])
  {
    for (auto const & name:
         config["additional_fields_p0"].as<std::vector<std::string>>())
    {
      auto const & [it, success] =
          fieldsP0_.emplace(name, proxpde::FEVar<FESpaceP0_T>{name, feSpaceP0_});
      assert(success);
      it->second.data = proxpde::Vec::Zero(feSpaceP0_.dof.size);
    }
  }

  // init bcs
  for (auto const & bc: config["bcs"])
  {
    auto const label = bc["label"].as<proxpde::marker_T>();
    auto const value = bc["value"].as<double>();
    bcs_.emplace_back(
        feSpace_, label, [value](proxpde::Vec3 const &) { return value; });
  }

  // init io
  if (config["output"])
  {
    outputPrefix_ = config["output"].as<std::string>();
  }
  io_.init(feSpace_, outputPrefix_ / "T");
  ioP0_.init(feSpaceP0_, outputPrefix_ / "fields");

  // coupling
  switch (equationType_)
  {
  case EQN_TYPE::HEAT:
  {
    couplingExport_.push_back("Tcfd");
    // couplingImport_.push_back("vel");
    break;
  }
  case EQN_TYPE::HEAT_COUPLED:
  {
    couplingExport_.push_back("T");
    couplingImport_.push_back("Tcfd");
    break;
  }
  case EQN_TYPE::HEAT_BUOYANT:
  {
    couplingExport_.push_back("T");
    couplingImport_.push_back("vel");
    break;
  }
  default:
    abort();
  }
}

uint ProblemProXPDEHeat::solve()
{
  auto const numIters = ProblemProXPDE::solve();

  // update local var
  T_.data = u_;

  return numIters;
}

void ProblemProXPDEHeat::print()
{
  if (it % printStep_ == 0)
  {
    // TODO: fix print of fields
    std::vector<proxpde::FEVar<FESpace_T>> fieldsP1 = {T_};
    for (auto & [_, value]: fieldsP1_)
    {
      fieldsP1.push_back(value);
    }
    io_.print(fieldsP1, time);
    std::vector<proxpde::FEVar<FESpaceP0_T>> fieldsP0;
    for (auto & [_, value]: fieldsP0_)
    {
      fieldsP0.push_back(value);
    }
    ioP0_.print(fieldsP0, time);
  }
}

namespace
{
template <typename FESpace>
std::vector<double> getData(proxpde::Vec const & u, FESpace const & feSpace)
{
  std::vector<double> data(u.size());
  for (auto const & e: feSpace.mesh->elementList)
    for (uint k = 0; k < FESpace::RefFE_T::numDOFs; ++k)
      for (uint d = 0; d < FESpace::dim; d++)
      {
        auto const id = FESpace::dim * e.pts[k]->id + d;
        data[id] = u[feSpace.dof.getId(e.id, k, d)];
      }
  return data;
}

template <typename FESpace>
void setData(FieldCoupling const & field, proxpde::Vec & u, FESpace const & feSpace)
{
  u.resize(field.size());
  for (auto const & e: feSpace.mesh->elementList)
    for (uint k = 0; k < FESpace::RefFE_T::numDOFs; ++k)
      for (uint d = 0; d < FESpace::dim; d++)
        u[feSpace.dof.getId(e.id, k, d)] = field[FESpace::dim * e.pts[k]->id + d];
}
} // namespace

std::unique_ptr<FieldCoupling> ProblemProXPDEHeat::initFieldCoupling(
    COUPLING_TYPE type, std::string_view name, MeshCoupling const * mesh)
{
  auto field = FieldCoupling::build(type);
  // check if the coupling field is the variable, T
  if (name == "T")
  {
    field->init(name, mesh, SUPPORT_TYPE::ON_NODES, NATURE_TYPE::INTENSIVE_MAXIMUM);
    auto const data = getData(T_.data, feSpace_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    field->initIO(outputPrefix_);
    return field;
  }
  // check if we are coupling the velocity
  if (name == "vel")
  {
    field->init(name, mesh, SUPPORT_TYPE::ON_NODES, NATURE_TYPE::INTENSIVE_MAXIMUM);
    auto const data = getData(vel_.data, feSpaceVel_);
    field->setValues({data.data(), data.data() + data.size()}, FESpaceVel_T::dim);
    field->initIO(outputPrefix_);
    return field;
  }
  // otherwise, it can be in fieldsP1_
  if (fieldsP1_.contains(std::string{name}))
  {
    field->init(name, mesh, SUPPORT_TYPE::ON_NODES, NATURE_TYPE::INTENSIVE_MAXIMUM);
    auto const data = getData(fieldsP1_.at(std::string{name}).data, feSpace_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    field->initIO(outputPrefix_);
    return field;
  }
  // otherwise, it can be in fieldsP0_
  if (fieldsP0_.contains(std::string{name}))
  {
    field->init(name, mesh, SUPPORT_TYPE::ON_CELLS, NATURE_TYPE::INTENSIVE_MAXIMUM);
    auto const data = getData(fieldsP0_.at(std::string{name}).data, feSpaceP0_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    field->initIO(outputPrefix_);
    return field;
  }

  fmt::println(stderr, "field {} not found", name);
  std::abort();
  return field;
}

void ProblemProXPDEHeat::setFieldData(FieldCoupling * field)
{
  // check if the coupling field is the variable, T
  if (field->name_ == "T")
  {
    auto const data = getData(T_.data, feSpace_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    return;
  }
  // check if we are coupling the velocity
  if (field->name_ == "vel")
  {
    auto const data = getData(vel_.data, feSpaceVel_);
    field->setValues({data.data(), data.data() + data.size()}, FESpaceVel_T::dim);
    return;
  }
  // otherwise, it can be in fieldsP1_
  if (fieldsP1_.contains(std::string{field->name_}))
  {
    auto const data = getData(fieldsP1_.at(std::string{field->name_}).data, feSpace_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    return;
  }
  // otherwise, it can be in fieldsP0_
  if (fieldsP0_.contains(std::string{field->name_}))
  {
    auto const data = getData(fieldsP0_.at(std::string{field->name_}).data, feSpaceP0_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    return;
  }

  fmt::println(stderr, "variable {} not found", field->name_);
  std::abort();
}

void ProblemProXPDEHeat::getFieldData(FieldCoupling const & field)
{
  // check if the coupling field is the variable, T
  if (field.name_ == "T")
  {
    setData(field, T_.data, feSpace_);
    return;
  }
  // check if the coupling field is the velocity
  if (field.name_ == "vel")
  {
    setData(field, vel_.data, feSpaceVel_);
    return;
  }
  // otherwise, it can be in fieldsP1_
  if (fieldsP1_.contains(std::string{field.name_}))
  {
    setData(field, fieldsP1_.at(std::string{field.name_}).data, feSpace_);
    return;
  }
  // otherwise, it can be in fieldsP0_
  if (fieldsP0_.contains(std::string{field.name_}))
  {
    setData(field, fieldsP0_.at(std::string{field.name_}).data, feSpaceP0_);
    return;
  }

  fmt::println(stderr, "variable {} not found", field.name_);
  std::abort();
}

// =====================================================================

struct AssemblyNS: public ProblemProXPDE::Assembly
{
  auto evaluate(ProblemProXPDE * pParent, proxpde::Builder<> & b) -> void override
  {
    auto p = dynamic_cast<ProblemProXPDENS *>(pParent);
    // assembly lhs
    b.buildLhs(
        std::tuple{
            proxpde::AssemblyMass{1. / p->dt_, p->feSpaceVel_},
            proxpde::AssemblyTensorStiffness{p->viscosity_, p->feSpaceVel_},
            proxpde::AssemblyAdvection{p->vel_, p->feSpaceVel_},
        },
        p->bcsVel_);
    b.buildCoupling(
        proxpde::AssemblyGrad{-1.0, p->feSpaceVel_, p->feSpaceP_},
        p->bcsVel_,
        p->bcsP_);
    b.buildCoupling(
        proxpde::AssemblyDiv{-1.0, p->feSpaceP_, p->feSpaceVel_}, p->bcsP_, p->bcsVel_);
    b.closeMatrix();

    // assembly rhs
    b.buildRhs(
        std::tuple{
            proxpde::AssemblyProjection{1. / p->dt_, p->vel_.data, p->feSpaceVel_}},
        p->bcsVel_);
  }
};

struct AssemblyNSBuoyant: public ProblemProXPDE::Assembly
{
  auto evaluate(ProblemProXPDE * pParent, proxpde::Builder<> & b) -> void override
  {
    auto p = dynamic_cast<ProblemProXPDENS *>(pParent);
    // update coupling
    proxpde::Vec T{p->feSpaceP_.dof.size};
    // p->getDataMED("T", T, p->feSpaceP_);

    // lhs terms
    auto const timeDer = proxpde::AssemblyMass{1. / p->dt_, p->feSpaceVel_};
    auto const advection = proxpde::AssemblyAdvection{p->vel_, p->feSpaceVel_};
    auto const diffusion =
        proxpde::AssemblyTensorStiffness{p->viscosity_, p->feSpaceVel_};
    auto const grad = proxpde::AssemblyGrad{-1.0, p->feSpaceVel_, p->feSpaceP_};
    auto const div = proxpde::AssemblyDiv{-1.0, p->feSpaceP_, p->feSpaceVel_};

    // rhs terms
    auto const timeDerOld =
        proxpde::AssemblyProjection{1. / p->dt_, p->vel_.data, p->feSpaceVel_};

    auto const beta = 1e+3;
    // auto const g = proxpde::Vec2{0.0, -1.0};
    auto const Tref = 0.0;
    using FESpaceBsq_T = proxpde::Vector_T<ProblemProXPDENS::FESpaceP_T, 2U>;
    auto feSpaceBsq = FESpaceBsq_T{p->mesh_};
    proxpde::Vec bsq = beta * (T.array() - Tref);
    // auto tmp2 = tmp * g.transpose();
    proxpde::Vec bsqVec = proxpde::Vec::Zero(bsq.size() * 2U);
    proxpde::setComponent(bsqVec, feSpaceBsq, bsq, p->feSpaceP_, 1u);
    auto const boussinesq =
        proxpde::AssemblyProjection{-1.0, bsqVec, feSpaceBsq, p->feSpaceVel_};

    // assembly lhs
    b.buildLhs(std::tuple{timeDer, advection, diffusion}, p->bcsVel_);
    b.buildCoupling(grad, p->bcsVel_, p->bcsP_);
    b.buildCoupling(div, p->bcsP_, p->bcsVel_);
    b.closeMatrix();

    // assembly rhs
    b.buildRhs(std::tuple{timeDerOld, boussinesq}, p->bcsVel_);
  }
};

ProblemProXPDENS::ProblemProXPDENS()
{
  [[maybe_unused]] auto const & [itNS, successNS] =
      assemblies_.emplace(EQN_TYPE::NS, new AssemblyNS);
  assert(successNS);

  [[maybe_unused]] auto const & [itNSB, successNSB] =
      assemblies_.emplace(EQN_TYPE::NS_BOUSSINESQ, new AssemblyNSBuoyant);
  assert(successNSB);
}

void ProblemProXPDENS::setup(Problem::ConfigList_T const & configs)
{
  // get configuration from file
  // TODO: string conversion is required for python bindings, manage other variant types
  // as errors
  auto const & configFileVariant = configs.at("config_file");
  auto const configFile =
      std::holds_alternative<std::filesystem::path>(configFileVariant)
          ? std::get<std::filesystem::path>(configFileVariant).string()
          : std::get<std::string>(configFileVariant);
  proxpde::ParameterDict config = YAML::LoadFile(configFile);

  name_ = config["name"].as<std::string>();

  // init mesh from configuration
  proxpde::readMesh(mesh_, config["mesh"]);

  // init equation
  equationType_ = str2eqn(config["equation_type"].as<std::string>());

  // init time related members
  time = 0.0;
  dt_ = config["dt"].as<double>();
  finalTime_ = config["final_time"].as<double>();

  // init fe related members
  feSpaceVel_.init(mesh_);
  feSpaceP_.init(mesh_, 2U * feSpaceVel_.dof.size);
  vel_.init("vel", feSpaceVel_);
  auto const velInitValue = config["vel_init_value"].as<proxpde::Vec2>();
  vel_ << velInitValue;

  p_.init("p", feSpaceP_);
  auto const pInitValue = config["p_init_value"].as<double>();
  p_.data = proxpde::Vec::Constant(feSpaceP_.dof.size, pInitValue);

  u_ = proxpde::Vec::Zero(2U * feSpaceVel_.dof.size + feSpaceP_.dof.size);
  u_.head(2U * feSpaceVel_.dof.size) = vel_.data;

  // init physical constants and source
  viscosity_ = config["viscosity"].as<double>();

  // init additional fields
  if (config["additional_fields_p1"])
  {
    for (auto const & name:
         config["additional_fields_p1"].as<std::vector<std::string>>())
    {
      auto const & [it, success] =
          fieldsP1_.emplace(name, proxpde::FEVar<FESpaceP_T>{name, feSpaceP_});
      assert(success);
      it->second.data = proxpde::Vec::Zero(feSpaceP_.dof.size);
    }
  }
  if (config["additional_fields_p0"])
  {
    for (auto const & name:
         config["additional_fields_p0"].as<std::vector<std::string>>())
    {
      auto const & [it, success] =
          fieldsP0_.emplace(name, proxpde::FEVar<FESpaceP0_T>{name, feSpaceP0_});
      assert(success);
      it->second.data = proxpde::Vec::Zero(feSpaceP0_.dof.size);
    }
  }

  // init bcs
  for (auto const & bc: config["bcs"])
  {
    auto const label = bc["label"].as<proxpde::marker_T>();
    auto const value = bc["value"].as<proxpde::Vec2>();
    bcsVel_.emplace_back(
        feSpaceVel_, label, [value](proxpde::Vec3 const &) { return value; });
  }

  // init io
  ioVel_.init(feSpaceVel_, "output_" + name_ + "/vel");
  ioP_.init(feSpaceP_, "output_" + name_ + "/p");

  // coupling
  feSpaceVelQ1_.init(mesh_);
  projectorQ2Q1_.init(feSpaceVelQ1_, feSpaceVel_);
  projectorQ2Q1_.setRhs(vel_.data);
  velQ1_.init("velQ1", feSpaceVelQ1_);
  velQ1_.data = projectorQ2Q1_.apply();
  couplingExport_.push_back("vel");

  switch (equationType_)
  {
  case EQN_TYPE::NS:
  {
    break;
  }
  case EQN_TYPE::NS_BOUSSINESQ:
  {
    // couplingImport_.push_back("T");
    // initFieldMED("T", "output_" + name_ + "/");
    // getField("T")->setValues(0.0, feSpaceP_.dof.size);
    break;
  }
  default:
    abort();
  }
}

uint ProblemProXPDENS::solve()
{
  auto const numIters = ProblemProXPDE::solve();

  // update local vars
  vel_.data = u_.head(2U * feSpaceVel_.dof.size);
  p_.data = u_.tail(feSpaceP_.dof.size);

  // update coupling vars
  projectorQ2Q1_.setRhs(vel_.data);
  velQ1_.data = projectorQ2Q1_.apply();

  return numIters;
}

void ProblemProXPDENS::print()
{
  if (it % printStep_ == 0)
  {
    ioVel_.print({vel_}, time);
    ioP_.print({p_}, time);
  }
}

std::unique_ptr<FieldCoupling> ProblemProXPDENS::initFieldCoupling(
    COUPLING_TYPE type, std::string_view name, MeshCoupling const * mesh)
{
  auto field = FieldCoupling::build(type);
  // check if the coupling field is the variable p
  if (name == "p")
  {
    field->init(name, mesh, SUPPORT_TYPE::ON_NODES, NATURE_TYPE::INTENSIVE_MAXIMUM);
    auto const data = getData(p_.data, feSpaceP_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    field->initIO(outputPrefix_);
    return field;
  }
  // check if the coupling field is the velocity
  if (name == "vel")
  {
    field->init(name, mesh, SUPPORT_TYPE::ON_NODES, NATURE_TYPE::INTENSIVE_MAXIMUM);
    auto const data = getData(velQ1_.data, feSpaceVelQ1_);
    field->setValues({data.data(), data.data() + data.size()}, FESpaceVel_T::dim);
    field->initIO(outputPrefix_);
    return field;
  }
  // otherwise, it can be in fieldsP1_
  if (fieldsP1_.contains(std::string{name}))
  {
    field->init(name, mesh, SUPPORT_TYPE::ON_NODES, NATURE_TYPE::INTENSIVE_MAXIMUM);
    auto const data = getData(fieldsP1_.at(std::string{name}).data, feSpaceP_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    field->initIO(outputPrefix_);
    return field;
  }
  // otherwise, it can be in fieldsP0_
  if (fieldsP0_.contains(std::string{name}))
  {
    field->init(name, mesh, SUPPORT_TYPE::ON_CELLS, NATURE_TYPE::INTENSIVE_MAXIMUM);
    auto const data = getData(fieldsP0_.at(std::string{name}).data, feSpaceP0_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    field->initIO(outputPrefix_);
    return field;
  }

  fmt::println(stderr, "field {} not found", name);
  std::abort();
  return field;
}

void ProblemProXPDENS::setFieldData(FieldCoupling * field)
{ // check if the coupling field is the variable p
  if (field->name_ == "p")
  {
    auto const data = getData(p_.data, feSpaceP_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    return;
  }
  // check if the coupling field is the velocity
  if (field->name_ == "vel")
  {
    auto const data = getData(velQ1_.data, feSpaceVelQ1_);
    field->setValues({data.data(), data.data() + data.size()}, FESpaceVel_T::dim);
    return;
  }
  // otherwise, it can be in fieldsP1_
  if (fieldsP1_.contains(std::string{field->name_}))
  {
    auto const data = getData(fieldsP1_.at(std::string{field->name_}).data, feSpaceP_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    return;
  }
  // otherwise, it can be in fieldsP0_
  if (fieldsP0_.contains(std::string{field->name_}))
  {
    auto const data = getData(fieldsP0_.at(std::string{field->name_}).data, feSpaceP0_);
    field->setValues({data.data(), data.data() + data.size()}, 1u);
    return;
  }

  fmt::println(stderr, "variable {} not found", field->name_);
  std::abort();
}
void ProblemProXPDENS::getFieldData(FieldCoupling const & field)
{ // check if the coupling field is the variable p
  if (field.name_ == "p")
  {
    setData(field, p_.data, feSpaceP_);
    return;
  }
  // check if the coupling field is the velocity
  if (field.name_ == "vel")
  {
    setData(field, vel_.data, feSpaceVel_);
    return;
  }
  // otherwise, it can be in fieldsP1_
  if (fieldsP1_.contains(std::string{field.name_}))
  {
    setData(field, fieldsP1_.at(std::string{field.name_}).data, feSpaceP_);
    return;
  }
  // otherwise, it can be in fieldsP0_
  if (fieldsP0_.contains(std::string{field.name_}))
  {
    setData(field, fieldsP0_.at(std::string{field.name_}).data, feSpaceP0_);
    return;
  }

  fmt::println(stderr, "variable {} not found", field.name_);
  std::abort();
}

} // namespace cocoa
