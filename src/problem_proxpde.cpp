#include "problem_proxpde.hpp"

// std
#include <string_view>

// proxpde
#include <proxpde/assembly.hpp>
#include <proxpde/builder.hpp>

// local
#include "enums.hpp"
#include "field_med.hpp"
#include "mesh_med.hpp"

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

template <typename Mesh>
void ProblemProXPDE::initMeshMED(std::string_view name, Mesh const & mesh)
{
  using Elem_T = Mesh::Elem_T;

  // coords format: x_0, y_0, z_0, x_1, ...
  std::vector<double> coords(mesh.pointList.size() * 3);
  for (auto const & p: mesh.pointList)
  {
    coords[3 * p.id] = p.coord[0];
    coords[3 * p.id + 1] = p.coord[1];
    coords[3 * p.id + 2] = p.coord[2];
  }

  // conn format: elem0_numpts, id_0, id_1, ..., elem1_numpts, ...
  std::vector<uint> conn(mesh.elementList.size() * (Elem_T::numPts + 1));
  for (auto const & elem: mesh.elementList)
  {
    conn[(Elem_T::numPts + 1) * elem.id] = MEDCellTypeToIKCell(MED_CELL_TYPE::QUAD4);
    for (uint p = 0; p < Elem_T::numPts; p++)
    {
      conn[(Elem_T::numPts + 1) * elem.id + p + 1] = elem.pts[p]->id;
    }
  }

  // offsets format: sum_0^k elemk_numpts + 1,
  std::vector<uint> offsets(mesh.elementList.size() + 1);
  offsets[0] = 0;
  for (uint e = 0; e < mesh.elementList.size(); e++)
  {
    offsets[e + 1] = offsets[e] + (Elem_T::numPts + 1);
  }

  meshCoupling_.reset(new MeshMED);
  meshCoupling_->init(name, Elem_T::dim, 3U, coords, conn, offsets);
}

void ProblemProXPDE::initFieldMED(std::string_view fieldName, std::string_view path)
{
  auto [kvPair, success] = fieldsCoupling_.emplace(fieldName, new FieldMED);
  assert(success);
  kvPair->second->init(fieldName, meshCoupling_.get());
  std::string filename = std::string{path} + "_med.";
  kvPair->second->initIO(filename);
}

template <typename FESpace>
void ProblemProXPDE::setDataMED(
    std::string_view fieldName, proxpde::Vec const & u, FESpace const & feSpace)
{
  std::vector<double> data(u.size());
  for (auto const & e: feSpace.mesh->elementList)
  {
    for (uint k = 0; k < FESpace::RefFE_T::numDOFs; ++k)
    {
      for (uint d = 0; d < FESpace::dim; d++)
      {
        data[FESpace::dim * e.pts[k]->id + d] = u[feSpace.dof.getId(e.id, k, d)];
      }
    }
  }
  getField(fieldName)->setValues(data, feSpace.dim);
}

template <typename FESpace>
void ProblemProXPDE::getDataMED(
    std::string_view fieldName, proxpde::Vec & u, FESpace const & feSpace)
{
  FieldCoupling * f = getField(fieldName);
  u.resize(f->size());

  for (auto const & e: feSpace.mesh->elementList)
  {
    for (uint k = 0; k < FESpace::RefFE_T::numDOFs; ++k)
    {
      for (uint d = 0; d < FESpace::dim; d++)
      {
        u[feSpace.dof.getId(e.id, k, d)] =
            *(f->dataPtr() + FESpace::dim * e.pts[k]->id + d);
      }
    }
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
}

bool ProblemProXPDE::run() { return time < finalTime_; }

void ProblemProXPDE::solve()
{
  fmt::print("{} - time = {:.6e}, dt = {:.6e}\n", name_, time, dt_);

  // update old values
  uOld_ = u_;

  // build matrix and rhs
  proxpde::Builder<> builder{this->size()};
  ProblemProXPDE::equationMap_.at(equationType_)(this, builder);

  // solve
  proxpde::LUSolver solver;
  solver.compute(builder.A);
  u_ = solver.solve(builder.b);
}

// =====================================================================

void ProblemProXPDEHeat::setup(ParamList_T const & params)
{
  // get configuration from file
  std::string const filename = params.at("config_file");
  proxpde::ParameterDict config = YAML::LoadFile(filename);

  // init mesh from configuration
  proxpde::readMesh(mesh_, config["mesh"]);
  name_ = config["name"].as<std::string>();
  initMeshMED(name_, mesh_);

  // init coupling
  // TODO: set from file when more coupling strategies are available
  couplingType_ = COUPLING_TYPE::MEDCOUPLING;

  // init equation
  equationType_ = str2eqn(config["equation_type"].as<std::string>());

  // init time related members
  time = 0.0;
  dt_ = config["dt"].as<double>();
  finalTime_ = config["final_time"].as<double>();

  // init fe related members
  feSpace_.init(mesh_);
  T_.init("T", feSpace_);
  double const initValue = config["init_value"].as<double>();
  T_ << initValue;
  u_ = T_.data;

  // init physical constants and operating conditions
  alpha_ = config["alpha"].as<double>();
  q_.init("source", feSpace_);
  q_ << config["q"].as<double>();
  feSpaceVel_.init(mesh_);
  vel_.init("vel", feSpaceVel_);
  vel_ << config["vel"].as<proxpde::Vec2>();

  // init bcs
  for (auto const & bc: config["bcs"])
  {
    auto const label = bc["label"].as<proxpde::marker_T>();
    auto const value = bc["value"].as<double>();
    bcs_.emplace_back(
        feSpace_, label, [value](proxpde::Vec3 const &) { return value; });
  }

  // init io
  io_.init(feSpace_, "output_" + name_ + "/T");

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
  initFieldMED(couplingExport_[0], io_.filePath->string());
  setDataMED(couplingExport_[0], T_.data, feSpace_);

  for (auto const & name: couplingImport_)
  {
    if (name != "vel")
    {
      initFieldMED(name, "output_" + name_ + "/i_" + name);
      setDataMED(name, T_.data, feSpace_);
    }
  }

  initFieldMED("vel", "output_" + name_ + "/i_vel");
  setDataMED("vel", vel_.data, feSpaceVel_);
}

void ProblemProXPDEHeat::solve()
{
  ProblemProXPDE::solve();

  // update local var
  T_.data = u_;

  // update coupling var
  setDataMED(couplingExport_[0], T_.data, feSpace_);
  setDataMED("vel", vel_.data, *vel_.feSpace);
}

void ProblemProXPDEHeat::print()
{
  // print med first to avoid iter update
  getField(couplingExport_[0])->printVTK(time, io_.iter);
  getField("vel")->printVTK(time, io_.iter);

  io_.print(std::tuple{T_, q_}, time);
}

// =====================================================================

void ProblemProXPDENS::setup(ParamList_T const & params)
{
  // get configuration from file
  std::string const filename = params.at("config_file");
  proxpde::ParameterDict config = YAML::LoadFile(filename);

  // init mesh from configuration
  proxpde::readMesh(mesh_, config["mesh"]);
  name_ = config["name"].as<std::string>();
  initMeshMED(name_, mesh_);

  // init coupling
  // TODO: set from file when more coupling strategies are available
  couplingType_ = COUPLING_TYPE::MEDCOUPLING;

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
  initFieldMED("vel", "output_" + name_ + "/vel");
  setDataMED("vel", velQ1_.data, feSpaceVelQ1_);

  switch (equationType_)
  {
  case EQN_TYPE::NS:
  {
    break;
  }
  case EQN_TYPE::NS_BOUSSINESQ:
  {
    couplingImport_.push_back("T");
    initFieldMED("T", "output_" + name_ + "/" + "i_T");
    getField("T")->setValues(0.0, feSpaceP_.dof.size);
    break;
  }
  default:
    abort();
  }
}

void ProblemProXPDENS::solve()
{
  ProblemProXPDE::solve();

  // update local vars
  vel_.data = u_.head(2U * feSpaceVel_.dof.size);
  p_.data = u_.tail(feSpaceP_.dof.size);

  // update coupling vars
  projectorQ2Q1_.setRhs(vel_.data);
  velQ1_.data = projectorQ2Q1_.apply();
  setDataMED("vel", velQ1_.data, feSpaceVelQ1_);
}

void ProblemProXPDENS::print()
{
  // print med first to avoid iter update
  for (auto const & name: couplingExport_)
  {
    getField(name)->printVTK(time, ioVel_.iter);
  }
  for (auto const & name: couplingImport_)
  {
    getField(name)->printVTK(time, ioVel_.iter);
  }

  ioVel_.print(std::tuple{vel_}, time);
  ioP_.print(std::tuple{p_}, time);
}

// =====================================================================

std::unordered_map<
    EQN_TYPE,
    std::function<void(ProblemProXPDE *, proxpde::Builder<> &)>>
    ProblemProXPDE::equationMap_ = {
        {EQN_TYPE::HEAT,
         [](ProblemProXPDE * problem, proxpde::Builder<> & b)
         {
           auto p = dynamic_cast<ProblemProXPDEHeat *>(problem);

           // update coupling
           // p->getDataMED("vel", p->vel_.data, *p->vel_.feSpace);

           // assembly lhs
           b.buildLhs(
               std::tuple{
                   proxpde::AssemblyMass{1. / p->dt_, p->feSpace_},
                   proxpde::AssemblyStiffness{p->alpha_, p->feSpace_},
                   proxpde::AssemblyAdvection{p->vel_, p->feSpace_},
               },
               p->bcs_);
           b.closeMatrix();

           // assembly rhs
           b.buildRhs(
               std::tuple{
                   proxpde::AssemblyProjection{1. / p->dt_, p->uOld_, p->feSpace_},
                   proxpde::AssemblyProjection{1., p->q_.data, p->feSpace_}},
               p->bcs_);
         }},
        {EQN_TYPE::HEAT_COUPLED,
         [](ProblemProXPDE * problem, proxpde::Builder<> & b)
         {
           auto p = dynamic_cast<ProblemProXPDEHeat *>(problem);

           // update coupling
           p->getDataMED("vel", p->vel_.data, *p->vel_.feSpace);

           // lhs terms
           auto timeDer = proxpde::AssemblyMass{1. / p->dt_, p->feSpace_};
           auto diffusion = proxpde::AssemblyStiffness{p->alpha_, p->feSpace_};
           auto advection = proxpde::AssemblyAdvection{p->vel_, p->feSpace_};

           proxpde::Vec mask(p->feSpace_.dof.size);
           p->getDataMED("mask", mask, p->feSpace_);
           auto const kAmpli = 20.;
           proxpde::FEVar kAmpliVec("kAmpli", p->feSpace_);
           kAmpliVec.data = kAmpli * mask;
           auto feedbackLhs = proxpde::AssemblyMassFE{1.0, kAmpliVec, p->feSpace_};

           // rhs terms
           auto timeDerOld =
               proxpde::AssemblyProjection{1. / p->dt_, p->uOld_, p->feSpace_};
           auto source = proxpde::AssemblyProjection{1., p->q_.data, p->feSpace_};

           proxpde::Vec Tcfd(p->feSpace_.dof.size);
           p->getDataMED("T", Tcfd, p->feSpace_);
           Tcfd.array() *= mask.array();
           auto feedbackRhs = proxpde::AssemblyProjection{kAmpli, Tcfd, p->feSpace_};

           // assembly
           b.buildLhs(std::tuple{timeDer, diffusion, advection, feedbackLhs}, p->bcs_);
           b.closeMatrix();
           b.buildRhs(std::tuple{timeDerOld, source, feedbackRhs}, p->bcs_);
         }},
        {EQN_TYPE::HEAT_BUOYANT,
         [](ProblemProXPDE * problem, proxpde::Builder<> & b)
         {
           auto p = dynamic_cast<ProblemProXPDEHeat *>(problem);

           // update coupling
           p->getDataMED("vel", p->vel_.data, *p->vel_.feSpace);
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
           auto diffusion = proxpde::AssemblyStiffness{p->alpha_, p->feSpace_};
           auto advection = proxpde::AssemblyAdvection{p->vel_, p->feSpace_};

           // rhs terms
           auto timeDerOld =
               proxpde::AssemblyProjection{1. / p->dt_, p->uOld_, p->feSpace_};
           // auto source = proxpde::AssemblyProjection{1., p->q_.data, p->feSpace_};

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
         }},
        {EQN_TYPE::NS,
         [](ProblemProXPDE * problem, proxpde::Builder<> & b)
         {
           auto p = dynamic_cast<ProblemProXPDENS *>(problem);

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
               proxpde::AssemblyDiv{-1.0, p->feSpaceP_, p->feSpaceVel_},
               p->bcsP_,
               p->bcsVel_);
           b.closeMatrix();

           // assembly rhs
           b.buildRhs(
               std::tuple{proxpde::AssemblyProjection{
                   1. / p->dt_, p->vel_.data, p->feSpaceVel_}},
               p->bcsVel_);
         }},
        {EQN_TYPE::NS_BOUSSINESQ,
         [](ProblemProXPDE * problem, proxpde::Builder<> & b)
         {
           auto p = dynamic_cast<ProblemProXPDENS *>(problem);

           // update coupling
           proxpde::Vec T{p->feSpaceP_.dof.size};
           p->getDataMED("T", T, p->feSpaceP_);

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
           proxpde::setComponent(bsqVec, feSpaceBsq, bsq, p->feSpaceP_, 1U);
           auto const boussinesq =
               proxpde::AssemblyProjection{-1.0, bsqVec, feSpaceBsq, p->feSpaceVel_};

           // assembly lhs
           b.buildLhs(std::tuple{timeDer, advection, diffusion}, p->bcsVel_);
           b.buildCoupling(grad, p->bcsVel_, p->bcsP_);
           b.buildCoupling(div, p->bcsP_, p->bcsVel_);
           b.closeMatrix();

           // assembly rhs
           b.buildRhs(std::tuple{timeDerOld, boussinesq}, p->bcsVel_);
         }},
};
