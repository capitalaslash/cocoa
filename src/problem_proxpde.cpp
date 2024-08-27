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

void ProblemProXPDE::setup(ParamList_T const & params)
{
  // get configuration from file
  std::string const filename = params.at("config_file");
  proxpde::ParameterDict config = YAML::LoadFile(filename);

  // init mesh from configuration
  proxpde::readMesh(mesh_, config["mesh"]);
  name_ = config["name"].as<std::string>();
  initMeshMED(name_);

  // init coupling
  // TODO: set from file when more coupling strategies are available
  couplingType_ = COUPLING_TYPE::MEDCOUPLING;

  // init equation
  equationType_ = str2proxpdeeqn(config["equation_type"].as<std::string>());

  // init time related members
  time = 0.0;
  dt_ = config["dt"].as<double>();
  finalTime_ = config["final_time"].as<double>();

  // init fe related members
  feSpace_.init(mesh_);
  auto const varName = config["var_name"].as<std::string>();
  u_.init(varName, feSpace_);
  double const initValue = config["init_value"].as<double>();
  u_.data = proxpde::Vec::Constant(feSpace_.dof.size, initValue);

  // init physical constants and source
  diff_ = config["diff"].as<double>();
  q_.init("source", feSpace_);
  q_.data = proxpde::Vec::Constant(feSpace_.dof.size, config["q"].as<double>());

  // init bcs
  for (auto const & bc: config["bcs"])
  {
    auto const label = bc["label"].as<proxpde::marker_T>();
    auto const value = bc["value"].as<double>();
    bcs_.emplace_back(
        feSpace_, label, [value](proxpde::Vec3 const &) { return value; });
  }

  // init io
  io_.init(feSpace_, "output_" + name_ + "/u");
  initFieldMED(varName);
}

void ProblemProXPDE::initMeshMED(std::string_view name)
{
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

  meshCoupling_.reset(new MeshMED);
  meshCoupling_->init(name, Elem_T::dim, 3U, coords, conn, offsets);
}

void ProblemProXPDE::initFieldMED(std::string_view name)
{
  auto [kvPair, success] = fieldsCoupling_.emplace(name, new FieldMED);
  assert(success);
  kvPair->second->init(name, meshCoupling_.get());
  setDataMED(u_.data, kvPair->first);
  std::string filename = std::string{io_.filePath.value()} + "_med.";
  kvPair->second->initIO(filename);
}

void ProblemProXPDE::setDataMED(proxpde::Vec const & u, std::string_view name)
{
  std::vector<double> data(u.size());
  for (auto const & e: mesh_.elementList)
  {
    for (uint k = 0; k < FESpace_T::RefFE_T::numDOFs; ++k)
    {
      data[e.pts[k]->id] = u[feSpace_.dof.getId(e.id, k)];
    }
  }
  getField(name)->setValues(data);
}

void ProblemProXPDE::getDataMED(proxpde::Vec & u, std::string_view name)
{
  FieldCoupling * f = getField(name);
  u.resize(f->size());

  for (auto const & e: mesh_.elementList)
  {
    for (uint k = 0; k < FESpace_T::RefFE_T::numDOFs; ++k)
    {
      u[feSpace_.dof.getId(e.id, k)] = *(f->dataPtr() + e.pts[k]->id);
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
  uOld_ = u_.data;

  // build matrix and rhs
  proxpde::Builder<> builder{feSpace_.dof.size};
  ProblemProXPDE::equationMap.at(equationType_)(this, builder);
  // builder.buildLhs(
  //     std::tuple{
  //         proxpde::AssemblyMass{1. / dt_, feSpace_},
  //         proxpde::AssemblyStiffness{diff_, feSpace_},
  //     },
  //     bcs_);
  // builder.closeMatrix();
  // builder.buildRhs(
  //     std::tuple{proxpde::AssemblyProjection{1. / dt_, uOld_, feSpace_}}, bcs_);

  // solve
  proxpde::LUSolver solver;
  solver.compute(builder.A);
  u_.data = solver.solve(builder.b);
}

void ProblemProXPDE::print()
{
  // print med first to avoid iter update
  setDataMED(u_.data, u_.name);
  getField(u_.name)->printVTK(time, io_.iter);

  io_.print(std::tuple{u_, q_}, time);
}

std::unordered_map<
    PROXPDEEQN_TYPE,
    std::function<void(ProblemProXPDE *, proxpde::Builder<> &)>>
    ProblemProXPDE::equationMap = {
        {PROXPDEEQN_TYPE::HEAT,
         [](ProblemProXPDE * p, proxpde::Builder<> & b)
         {
           b.buildLhs(
               std::tuple{
                   proxpde::AssemblyMass{1. / p->dt_, p->feSpace_},
                   proxpde::AssemblyStiffness{p->diff_, p->feSpace_},
               },
               p->bcs_);
           b.closeMatrix();
           b.buildRhs(
               std::tuple{
                   proxpde::AssemblyProjection{1. / p->dt_, p->uOld_, p->feSpace_},
                   proxpde::AssemblyProjection{1., p->q_.data, p->feSpace_}},
               p->bcs_);
         }},
        {PROXPDEEQN_TYPE::HEAT_COUPLED,
         [](ProblemProXPDE * p, proxpde::Builder<> & b)
         {
           auto timeDer = proxpde::AssemblyMass{1. / p->dt_, p->feSpace_};
           auto diffusion = proxpde::AssemblyStiffness{p->diff_, p->feSpace_};
           auto const kAmpli = 20.;

           proxpde::Vec mask(p->feSpace_.dof.size);
           p->getDataMED(mask, "mask");
           proxpde::FEVar kAmpliVec("kAmpli", p->feSpace_);
           kAmpliVec.data = kAmpli * mask;
           auto feedbackLhs = proxpde::AssemblyMassFE{1.0, kAmpliVec, p->feSpace_};

           auto timeDerOld =
               proxpde::AssemblyProjection{1. / p->dt_, p->uOld_, p->feSpace_};
           auto source = proxpde::AssemblyProjection{1., p->q_.data, p->feSpace_};
           proxpde::Vec uExt(p->feSpace_.dof.size);
           p->getDataMED(uExt, "Tcfd");
           uExt.array() *= mask.array();
           auto feedbackRhs = proxpde::AssemblyProjection{kAmpli, uExt, p->feSpace_};

           b.buildLhs(std::tuple{timeDer, diffusion, feedbackLhs}, p->bcs_);
           b.closeMatrix();
           b.buildRhs(std::tuple{timeDerOld, source, feedbackRhs}, p->bcs_);
         }},
};
