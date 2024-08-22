#include "problem_proxpde.hpp"

// std
#include <string_view>

// proxpde
#include <proxpde/assembly.hpp>
#include <proxpde/builder.hpp>

// local
#include "field_med.hpp"
#include "mesh_med.hpp"

void ProblemProXPDE::setup(ParamList_T const & params)
{
  // get configuration from file
  std::string const filename = params.at("config_file");
  proxpde::ParameterDict config = YAML::LoadFile(filename);

  // init mesh from configuration
  proxpde::readMesh(mesh_, config["mesh"]);
  auto const name = config["name"].as<std::string>();
  initMeshMED(name);

  // init time related members
  time = 0.0;
  dt_ = config["dt"].as<double>();
  finalTime_ = config["final_time"].as<double>();

  // init fe related members
  feSpace_.init(mesh_);
  diff_ = config["diff"].as<double>();
  auto const varName = config["var_name"].as<std::string>();
  u_.init(varName, feSpace_);
  double const initValue = config["init_value"].as<double>();
  u_.data = proxpde::Vec::Constant(feSpace_.dof.size, initValue);

  // init bcs
  for (auto const & bc: config["bcs"])
  {
    auto const label = bc["label"].as<proxpde::marker_T>();
    auto const value = bc["value"].as<double>();
    bcs_.emplace_back(
        feSpace_, label, [value](proxpde::Vec3 const &) { return value; });
  }

  // init io
  io_.init(feSpace_, "output_" + name + "/u");
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

  meshCoupling_->init(name, Elem_T::dim, 3U, coords, conn, offsets);
}

void ProblemProXPDE::initFieldMED(std::string_view name)
{
  auto [kvPair, success] = fieldsCoupling_.emplace(u_.name, new FieldMED);
  assert(success);
  kvPair->second->init(name, meshCoupling_.get());
  updateFieldMED(kvPair->first);
  std::string filename = std::string{io_.filePath.value()} + "_med.";
  kvPair->second->initIO(filename);
}

void ProblemProXPDE::updateFieldMED(std::string_view name)
{
  std::vector<double> data(u_.data.size());
  for (auto const & e: mesh_.elementList)
  {
    for (uint k = 0; k < FESpace_T::RefFE_T::numDOFs; ++k)
    {
      data[e.pts[k]->id] = u_.data[feSpace_.dof.getId(e.id, k)];
    }
  }
  getField(name)->setValues(data);
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
  fmt::print("proxpde::solve(), time = {:.6e}, dt = {:.6e}\n", time, dt_);

  // update old values
  uOld_ = u_.data;

  // build matrix and rhs
  proxpde::Builder<> builder{feSpace_.dof.size};
  builder.buildLhs(
      std::tuple{
          proxpde::AssemblyMass{1. / dt_, feSpace_},
          proxpde::AssemblyStiffness{diff_, feSpace_},
      },
      bcs_);
  builder.closeMatrix();
  builder.buildRhs(
      std::tuple{proxpde::AssemblyProjection{1. / dt_, uOld_, feSpace_}}, bcs_);

  // solve
  proxpde::LUSolver solver;
  solver.compute(builder.A);
  u_.data = solver.solve(builder.b);
}

void ProblemProXPDE::print()
{
  // print med before to avoid iter update
  updateFieldMED(u_.name);
  getField(u_.name)->printVTK(time, io_.iter);

  io_.print(std::tuple{u_}, time);
}