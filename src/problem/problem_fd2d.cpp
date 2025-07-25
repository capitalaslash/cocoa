#include "problem/problem_fd2d.hpp"

// std
#include <cassert>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <unordered_map>

// libfmt
#include <fmt/core.h>
#include <fmt/ranges.h>

// local
#include "coupling/field_coupling.hpp"
#include "coupling/mesh_coupling.hpp"
#include "enums.hpp"
#include "problem/fdutils.hpp"

namespace cocoa
{

ProblemFD2D::ProblemFD2D(): Problem{PROBLEM_TYPE::FD2D}
{
  // register default assemblies
  assemblies_.emplace(EQN_TYPE::HEAT, [](ProblemFD2D * p) { p->assemblyHeat(); });
  assemblies_.emplace(
      EQN_TYPE::HEAT_COUPLED, [](ProblemFD2D * p) { p->assemblyHeatCoupled(); });
  assemblies_.emplace(EQN_TYPE::HEAT_OC, [](ProblemFD2D * p) { p->assemblyHeatOC(); });
}

ProblemFD2D::~ProblemFD2D()
{
  // erase possibly added assembly
  assemblies_.erase(EQN_TYPE::CUSTOM);
}

constexpr FD_BC_SIDE side2D(uint s)
{
  switch (s)
  {
  case 0u:
    return FD_BC_SIDE::LEFT;
  case 1u:
    return FD_BC_SIDE::RIGHT;
  case 2u:
    return FD_BC_SIDE::BOTTOM;
  case 3u:
    return FD_BC_SIDE::TOP;
  default:
    std::abort();
  }
}

void ProblemFD2D::setup(Problem::ConfigList_T const & configs)
{
  // default values
  name_ = "empty";
  // mesh
  MeshFD2D::Real_T start = {0.0, 0.0};
  MeshFD2D::Real_T end = {1.0, 1.0};
  MeshFD2D::Int_T nElems = {10u, 10u};
  // fields
  nVars_ = 1u;
  varNames_ = {"u"};
  std::vector<double> uInit(nVars_, 0.0);
  std::vector<double> qValue(nVars_, 1.0);
  // bcs
  bcs_.resize(nVars_);
  std::vector<std::array<FD_BC_TYPE, 4u>> bcTypes(nVars_);
  std::vector<std::array<double, 4u>> bcValues(nVars_);
  // time
  time = 0.0;
  finalTime_ = 1.0;
  dt_ = 0.1;
  // linear algebra
  maxIters_ = 1000u;
  tol_ = 1.e-6;
  uint nnz = 5u;

  // read configuration from file
  // TODO: string conversion is required for python bindings, manage other variant types
  // as errors
  auto const & configFileVariant = configs.at("config_file");
  auto const configFile =
      std::holds_alternative<std::filesystem::path>(configFileVariant)
          ? std::get<std::filesystem::path>(configs.at("config_file"))
          : std::filesystem::path{std::get<std::string>(configs.at("config_file"))};
  std::ifstream in(configFile, std::ios::in);
  if (!in)
  {
    fmt::println(stderr, "configuration file {} not found!", configFile.string());
    std::abort();
  }
  std::string buffer;
  while (std::getline(in, buffer, '\n'))
  {
    std::istringstream bufferStream{buffer};
    std::string token;
    while (std::getline(bufferStream, token, ' '))
    {
      if (token[0] == '#')
      {
        // this is a comment, ignore line
      }
      else if (token == "name:")
        bufferStream >> name_;
      else if (token == "debug:")
        bufferStream >> debug_;
      else if (token == "compute_cfl:")
        bufferStream >> computeCFL_;
      // mesh
      else if (token == "start:")
      {
        bufferStream >> start[0];
        bufferStream >> start[1];
      }
      else if (token == "end:")
      {
        bufferStream >> end[0];
        bufferStream >> end[1];
      }
      else if (token == "n_elems:")
      {
        bufferStream >> nElems[0];
        bufferStream >> nElems[1];
      }
      // fields
      else if (token == "n_vars:")
      {
        bufferStream >> nVars_;
        varNames_.resize(nVars_);
        uInit.resize(nVars_);
        qValue.resize(nVars_);
        bcs_.resize(nVars_);
        bcTypes.resize(nVars_);
        bcValues.resize(nVars_);
      }
      else if (token == "var_names:")
        for (uint v = 0u; v < nVars_; v++)
          bufferStream >> varNames_[v];
      else if (token == "initial_value:")
        for (uint v = 0u; v < nVars_; v++)
          bufferStream >> uInit[v];
      else if (token == "q:")
        for (uint v = 0u; v < nVars_; v++)
          bufferStream >> qValue[v];
      // params
      else if (token == "params:")
      {
        while (bufferStream)
        {
          std::string name = "";
          bufferStream >> name;
          if (name == "")
            break;
          bufferStream >> token;
          FD_PARAM_TYPE type = str2FDParamType(token);
          switch (type)
          {
          case FD_PARAM_TYPE::INTEGER:
          {
            uint value;
            bufferStream >> value;
            params_.set(name, value);
            break;
          }
          case FD_PARAM_TYPE::SCALAR:
          {
            double value;
            bufferStream >> value;
            params_.set(name, value);
            break;
          }
          case FD_PARAM_TYPE::VECTOR:
          {
            uint size;
            bufferStream >> size;
            assert(size > 0u);
            std::vector<double> value(size);
            for (uint k = 0u; k < size; k++)
              bufferStream >> value[k];
            params_.set(name, value);
            break;
          }
          default:
            fmt::println(stderr, "param type for {} not recognized", name);
            std::abort();
          }
        }
      }
      // time
      else if (token == "start_time:")
        bufferStream >> time;
      else if (token == "final_time:")
        bufferStream >> finalTime_;
      else if (token == "dt:")
        bufferStream >> dt_;
      // assembly
      else if (token == "assembly_name:")
      {
        bufferStream >> token;
        eqnType_ = str2eqn(token);
      }
      // la
      else if (token == "solver_type:")
      {
        bufferStream >> token;
        solverType_ = str2fdsolver(token);
      }
      else if (token == "max_iters:")
        bufferStream >> maxIters_;
      else if (token == "tol:")
        bufferStream >> tol_;
      else if (token == "nnz:")
        bufferStream >> nnz;
      // bcs
      else if (token == "bc_bottom:")
      {
        for (uint v = 0u; v < nVars_; v++)
        {
          bufferStream >> token;
          bcTypes[v][2u] = str2FDBCType(token);
          bufferStream >> bcValues[v][2u];
        }
      }
      else if (token == "bc_right:")
      {
        for (uint v = 0u; v < nVars_; v++)
        {
          bufferStream >> token;
          bcTypes[v][1u] = str2FDBCType(token);
          bufferStream >> bcValues[v][1u];
        }
      }
      else if (token == "bc_top:")
      {
        for (uint v = 0u; v < nVars_; v++)
        {
          bufferStream >> token;
          bcTypes[v][3u] = str2FDBCType(token);
          bufferStream >> bcValues[v][3u];
        }
      }
      else if (token == "bc_left:")
      {
        for (uint v = 0u; v < nVars_; v++)
        {
          bufferStream >> token;
          bcTypes[v][0u] = str2FDBCType(token);
          bufferStream >> bcValues[v][0u];
        }
      }
      // io
      else if (token == "print_step:")
        bufferStream >> printStep_;
      else if (token == "output_prefix:")
        bufferStream >> outputPrefix_;
      else if (token == "clean_output:")
        bufferStream >> cleanOutput_;
      else
      {
        fmt::println(stderr, "key {} invalid", token);
        bufferStream >> token;
      }
    }
  }
  fmt::println("{} - equation type: {}", name_, eqn2str(eqnType_));
  assert(eqnType_ == EQN_TYPE::NONE || assemblies_.contains(eqnType_));

  fmt::println("parameters: {}", params_);

  // mesh
  mesh_.init(start, end, nElems);

  // fields
  u_.resize(mesh_.nPts() * nVars_);
  uOld_.resize(mesh_.nPts() * nVars_);
  q_.resize(mesh_.nPts() * nVars_);
  std::vector<double> cValues(nVars_ * 2u, 0.0);
  if (params_.data_.contains("c"))
  {
    cValues = params_.get<FD_PARAM_TYPE::VECTOR>("c");
    assert(cValues.size() == nVars_ * 2u);
    if (nVars_ == 1u)
    {
      fields_.emplace("cx", VectorFD{mesh_.nPts(), cValues[0]});
      fields_.emplace("cy", VectorFD{mesh_.nPts(), cValues[1]});
    }
    else
    {
      for (uint v = 0u; v < nVars_; v++)
      {
        fields_.emplace(
            "cx" + std::to_string(v), VectorFD{mesh_.nPts(), cValues[0 + 2 * v]});
        fields_.emplace(
            "cy" + std::to_string(v), VectorFD{mesh_.nPts(), cValues[1 + 2 * v]});
      }
    }
  }
  for (uint v = 0u; v < nVars_; v++)
  {
    u_.setRange(0u + v * mesh_.nPts(), mesh_.nPts() + v * mesh_.nPts(), uInit[v]);
    uOld_.setRange(0u + v * mesh_.nPts(), mesh_.nPts() + v * mesh_.nPts(), uInit[v]);
    q_.setRange(0u + v * mesh_.nPts(), mesh_.nPts() + v * mesh_.nPts(), qValue[v]);

    for (uint s = 0u; s < 4u; s++)
      bcs_[v].data_[s] =
          FDBC(side2D(s), bcTypes[v][s], bcValues[v][s], mesh_.n_[1 - s / 2]);
  }

  // linear algebra
  m_.init(mesh_.nPts() * nVars_, nnz);
  rhs_.resize(mesh_.nPts() * nVars_);

  // io
  initOutput();
}

std::unique_ptr<MeshCoupling> ProblemFD2D::initMeshCoupling(
    COUPLING_TYPE type, COUPLING_SCOPE scope, Marker marker, std::string_view bdName)
{
  if (scope == COUPLING_SCOPE::VOLUME)
  {
    // coords format: x_0, y_0, z_0, x_1, ...
    std::vector<double> coords(mesh_.nPts() * 3);
    for (uint j = 0; j < mesh_.n_[1]; j++)
      for (uint i = 0; i < mesh_.n_[0]; i++)
      {
        uint const id = j * mesh_.n_[0] + i;
        auto const pt = mesh_.pt({i, j});
        coords[3 * id + 0U] = pt[0];
        coords[3 * id + 1u] = pt[1];
        coords[3 * id + 2U] = 0.0;
      }

    // conn format: elem0_numpts, id_0, id_1, ..., elem1_numpts, ...
    auto const nElems = mesh_.nElems();
    std::vector<uint> conn(nElems * (1 + 4));
    auto elemCount = 0U;
    for (uint j = 0; j < mesh_.n_[1] - 1; j++)
      for (uint i = 0; i < mesh_.n_[0] - 1; i++)
      {
        uint const id = j * mesh_.n_[0] + i;
        conn[(1 + 4) * elemCount] = 4; // MEDCellTypeToIKCell(MED_CELL_TYPE::QUAD4);
        conn[(1 + 4) * elemCount + 1] = id;
        conn[(1 + 4) * elemCount + 2] = id + 1;
        conn[(1 + 4) * elemCount + 3] = id + mesh_.n_[0] + 1;
        conn[(1 + 4) * elemCount + 4] = id + mesh_.n_[0];
        elemCount++;
      }

    // offsets format: sum_0^k elemk_numpts + 1,
    std::vector<uint> offsets(nElems + 1);
    offsets[0] = 0;
    for (uint k = 0; k < nElems; k++)
    {
      offsets[k + 1] = offsets[k] + (1 + 4);
    }

    auto meshCoupling = MeshCoupling::build(type);
    meshCoupling->init(
        "mesh_fd2d",
        COUPLING_SCOPE::VOLUME,
        markerNotSet,
        "",
        2u,
        2u,
        coords,
        conn,
        offsets);
    // meshCoupling->printVTK(outputPrefix_ / "mesh_fd2d");
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

std::unique_ptr<FieldCoupling> ProblemFD2D::initFieldCoupling(
    COUPLING_TYPE type, std::string_view name, MeshCoupling const * mesh)
{
  auto field = FieldCoupling::build(type);
  for (uint v = 0u; v < nVars_; v++)
  {
    if (varNames_[v] == name)
    {
      field->init(name, mesh, SUPPORT_TYPE::ON_NODES, NATURE_TYPE::INTENSIVE_MAXIMUM);
      field->setValues(
          {u_.data() + v * mesh_.nPts(), u_.data() + (v + 1) * mesh_.nPts()}, 1u);
      field->initIO(outputPrefix_);
      return field;
    }
  }
  // check if the coupling field is an additional field
  if (fields_.contains(std::string{name}))
  {
    field->init(name, mesh, SUPPORT_TYPE::ON_NODES, NATURE_TYPE::INTENSIVE_MAXIMUM);
    field->setValues(fields_.at(std::string{name}).data_, 1u);
    return field;
  }

  fmt::println(stderr, "coupling field {} not found", name);
  std::abort();
  return field;
}

void ProblemFD2D::setFieldData(FieldCoupling * field)
{
  // check if the field maps a variable
  for (auto v = 0u; v < nVars_; v++)
  {
    if (varNames_[v] == field->name_)
    {
      auto const start = u_.data() + v * mesh_.nPts();
      field->setValues({start, start + mesh_.nPts()});
      return;
    }
  }
  // otherwise, it must be stored in the additional fields
  field->setValues(fields_.at(field->name_).data_);
}

void ProblemFD2D::getFieldData(FieldCoupling const & field)
{
  // check if the field maps a variable
  for (auto v = 0u; v < nVars_; v++)
  {
    if (varNames_[v] == field.name_)
    {
      auto const start = u_.data() + v * mesh_.nPts();
      std::copy(field.dataPtr(), field.dataPtr() + field.size(), start);
      return;
    }
  }
  // otherwise, it must be stored in the additional fields
  std::copy(
      field.dataPtr(), field.dataPtr() + field.size(), fields_.at(field.name_).data());
}

void ProblemFD2D::initOutput()
{
  if (cleanOutput_ && std::filesystem::exists(outputPrefix_))
    for (const auto & entry: std::filesystem::directory_iterator(outputPrefix_))
      std::filesystem::remove_all(entry.path());
  std::filesystem::create_directories(outputPrefix_);
}

bool ProblemFD2D::run() const { return time < finalTime_; }

void ProblemFD2D::advance()
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

// TODO: static constexpr std::vector<uint> requires gcc >= 12
const std::vector<uint> sideDOF(std::array<uint, 2U> const & n, FD_BC_SIDE const side)
{
  std::vector<uint> dofList;

  switch (side)
  {
  case FD_BC_SIDE::LEFT:
  {
    dofList.resize(n[1]);
    for (uint k = 0u; k < dofList.size(); k++)
      dofList[k] = k * n[0];
    break;
  }
  case FD_BC_SIDE::RIGHT:
  {
    dofList.resize(n[1]);
    for (uint k = 0u; k < dofList.size(); k++)
      dofList[k] = (k + 1) * n[0] - 1;
    break;
  }
  case FD_BC_SIDE::BOTTOM:
  {
    dofList.resize(n[0]);
    std::iota(dofList.begin(), dofList.end(), 0u);
    break;
  }
  case FD_BC_SIDE::TOP:
  {
    dofList.resize(n[0]);
    std::iota(dofList.begin(), dofList.end(), n[0] * (n[1] - 1));
    break;
  }
  default:
    std::abort();
  }

  return dofList;
}

// static constexpr std::array<uint, 4U> cornerDOF(std::array<uint, 2U> const & n)
// {
//   return std::array<uint, 4U>{{
//       0U,                // bottom-left
//       n[0] - 1,          // bottom-right
//       n[0] * n[1] - 1,   // top-right
//       n[0] * (n[1] - 1), // top-left
//   }};
// }

// static constexpr std::array<std::array<uint, 2U>, 4U> cornerSides = {{
//     {{0, 3}}, // bottom-left
//     {{1, 0}}, // bottom-right
//     {{2, 1}}, // top-right
//     {{3, 2}}, // top-left
// }};

// static constexpr int cornerOffset(std::array<uint, 2U> const & n, uint k)
// {
//   switch (k)
//   {
//   case 0U: // bottom-left
//     return n[0] + 1;
//   case 1u: // bottom-right
//     return n[0] - 1;
//   case 2U: // top-right
//     return -n[0] - 1;
//   case 3U: // top-left
//     return -n[0] + 1;
//   default:
//     std::abort();
//   }
//   return 0;
// }

// static constexpr std::pair<uint, uint> cornerEnd(std::array<uint, 2U> const & n,
// uint k)
// {
//   switch (k)
//   {
//   case 0U:
//     return {0, 0}; // bottom left
//   case 1u:
//     return {0, n[0] - 1}; // bottom right
//   case 2U:
//     return {n[0] - 1, n[1] - 1}; // top left
//   case 3U:
//     return {n[1] - 1, 0}; // top right
//   default:
//     std::abort();
//   }
//   return {0, 0};
// }

uint ProblemFD2D::solve()
{
  fmt::println("\n===");
  fmt::println("{}, time = {:.6e}, dt = {:.6e}", name_, time, dt_);

  // TODO: improve CFL evaluation by using better estimation of cell diameter
  if (computeCFL_)
  {
    auto const cx = fields_.at("cx");
    auto const cy = fields_.at("cy");
    double maxCFL = 0.0;
    for (uint k = 0U; k < u_.size(); k++)
    {
      double const cLocal = std::sqrt(cx[k] * cx[k] + cy[k] * cy[k]);
      maxCFL = std::max(cLocal * dt_ / std::min(mesh_.h_[0], mesh_.h_[1]), maxCFL);
    }
    fmt::println("maxCFL: {:.6e}", maxCFL);
  }
  // update
  uOld_ = u_;

  // assembly
  assemblies_.at(eqnType_)(this);

  // fmt::println("bc b: {}", bcs_[0].bottom());
  // fmt::println("bc r: {}", bcs_[0].right());
  // fmt::println("bc t: {}", bcs_[0].top());
  // fmt::println("bc l: {}", bcs_[0].left());

  std::array<double, 4U> const hSide = {
      mesh_.h_[1], mesh_.h_[1], mesh_.h_[0], mesh_.h_[0]};
  for (uint v = 0u; v < nVars_; v++)
  {
    // bc: Neumann sides
    for (uint s = 0u; s < 4u; s++)
    {
      auto const & bc = bcs_[v].data_[s];
      auto const dofList = sideDOF(mesh_.n_, side2D(s));

      if (bc.type == FD_BC_TYPE::NEUMANN)
      {
        for (uint k = 0U; k < dofList.size(); k++)
        {
          uint const dof = dofList[k] + v * mesh_.nPts();
          // sign: incoming flux is positive
          // (u_in - u_out) / 2h = A
          // u_out = u_in - 2 h A
          // u_in part implemented in assembly
          rhs_.add(dof, 2.0 * hSide[s] * bc.values[k] * bc.ghostValues[k]);
        }
      }
    }

    // bc: Dirichlet sides
    for (uint s = 0U; s < 4U; s++)
    {
      auto const & bc = bcs_[v].data_[s];

      if (bc.type == FD_BC_TYPE::DIRICHLET)
      {
        auto const dofList = sideDOF(mesh_.n_, side2D(s));
        for (uint k = 0U; k < dofList.size(); k++)
        {
          uint const dof = dofList[k] + v * mesh_.nPts();
          m_.clearRow(dof);
          m_.add(dof, dof, 1.0);
          rhs_.set(dof, bc.values[k]);
        }
      }
    }

    // TODO: manage Neumann/Neumann corners?
    // for (uint k = 0; k < 4U; k++)
    // {
    //   auto const dof = cornerDOF(n_)[k];
    //   m_.add(dof, dof, 1.0);
    //   if (bcs_[cornerSides[k][0]].type == FD_BC_TYPE::NEUMANN &&
    //   bcs_[cornerSides[k][1]].type == FD_BC_TYPE::NEUMANN)
    //   {
    //     m_.add(dof, dof + cornerOffset(n_, k), -1.0);
    //     rhs_[dof] =
    //         bcs_[cornerSides[k][0]].values[] * h_[1] +
    //         bcs_[cornerSides[k][1]].values[]
    //         * h_[0];
    //   }
    // }
  }

  m_.close();

  // pre-solve
  if (preSolveFun_)
    (*preSolveFun_)(this);

  if (debug_)
    m_.print_sparsity_pattern("fd2d_mat.dat");

  // solve
  auto const [numIters, residual] =
      solvers_.at(solverType_)(m_, rhs_, u_, tol_, maxIters_);
  fmt::print("num iters: {:4d}, ", numIters);
  fmt::println("relative residual: {:.8e}", residual);

  if (debug_)
  {
    fmt::println("matrix: {}", m_);
    fmt::println("rhs: {}", rhs_);
    fmt::println("sol: {}", u_);
  }
  // clean up
  m_.clear();
  rhs_.zero();

  return numIters;
}

void ProblemFD2D::assemblyHeat()
{
  auto const h = mesh_.h_;
  auto const alpha = params_.get<FD_PARAM_TYPE::SCALAR>("alpha");
  auto const & cx = fields_.at("cx");
  auto const & cy = fields_.at("cy");

  for (uint j = 0u; j < mesh_.n_[1]; j++)
    for (uint i = 0u; i < mesh_.n_[0]; i++)
    {
      auto const id = i + j * mesh_.n_[0];

      // eqn:
      // du/dt + c du/dx - alpha * d^2u/dx^2 = q
      // discretization:
      // u_m / dt + cx * (u_r - u_l) / (2 * hx) + cy * (u_t - u_b) / (2 * hy)
      // - alpha (u_l - 2 * u_m + u_r) / hx^2 - alpha (u_b - 2 * u_m + u_t) / hy^2 =
      // uold_m / dt + q_m
      // grouping:
      // (1 / dt + 2 * alpha (1 / hx^2 + 1 / hy^2) * u_m
      // - (alpha / hx^2 - cx / hx) * u_l
      // - (alpha / hx^2 + cx / hx) * u_r
      // - (alpha / hy^2 - cy / hy) * u_b
      // - (alpha / hy^2 + cy / hy) * u_t
      // = u_m / dt + q_m

      // middle
      double const value =
          1. / dt_                                                 // time derivative
          + 2. * alpha * (1. / (h[0] * h[0]) + 1. / (h[1] * h[1])) // diffusion
          ;
      m_.add(id, id, value);

      // bottom
      double const valueBottom = -alpha / (h[1] * h[1]) // diffusion
                                 - 0.5 * cy[id] / h[1]  // advection
          ;
      if (j > 0u)
        m_.add(id, id - mesh_.n_[0], valueBottom);
      else
      {
        m_.add(id, id + mesh_.n_[0], valueBottom);
        bcs_[0].bottom().ghostValues.set(i, valueBottom);
      }

      // right
      double const valueRight = -alpha / (h[0] * h[0]) // diffusion
                                + 0.5 * cx[id] / h[0]  // advection
          ;
      if (i < mesh_.n_[0] - 1)
        m_.add(id, id + 1, valueRight);
      else
      {
        m_.add(id, id - 1, valueRight);
        bcs_[0].right().ghostValues.set(j, valueRight);
      }

      // top
      double const valueTop = -alpha / (h[1] * h[1]) // diffusion
                              + 0.5 * cy[id] / h[1]  // advection
          ;
      if (j < mesh_.n_[1] - 1)
        m_.add(id, id + mesh_.n_[0], valueTop);
      else
      {
        m_.add(id, id - mesh_.n_[0], valueTop);
        bcs_[0].top().ghostValues.set(i, valueTop);
      }

      // left
      double const valueLeft = -alpha / (h[0] * h[0]) // diffusion
                               - 0.5 * cx[id] / h[0]  // advection
          ;
      if (i > 0u)
        m_.add(id, id - 1, valueLeft);
      else
      {
        m_.add(id, id + 1, valueLeft);
        bcs_[0].left().ghostValues.set(j, valueLeft);
      }

      // rhs
      rhs_.set(
          id,
          uOld_[id] / dt_ // time derivative
              + q_[id]    // source
      );
    }
  m_.close();
}

void ProblemFD2D::assemblyHeatOC()
{
  auto const h = mesh_.h_;
  auto const eps = 0.5 * std::min(h[0], h[1]);
  auto const alpha = params_.get<FD_PARAM_TYPE::SCALAR>("alpha");
  auto const beta = params_.get<FD_PARAM_TYPE::SCALAR>("beta");
  auto const target = params_.get<FD_PARAM_TYPE::VECTOR>("target");
  auto const control = params_.get<FD_PARAM_TYPE::VECTOR>("control");

  for (uint j = 0u; j < mesh_.n_[1]; j++)
    for (uint i = 0u; i < mesh_.n_[0]; i++)
    {
      auto const idF = i + j * mesh_.n_[0];
      auto const idA = idF + mesh_.nPts();
      auto const pt = mesh_.pt({i, j});

      // forward
      // middle
      double const valueF =
          1. / dt_                                                 // time derivative
          + 2. * alpha * (1. / (h[0] * h[0]) + 1. / (h[1] * h[1])) // diffusion
          ;
      m_.add(idF, idF, valueF);

      // bottom
      double const valueBottomF = -alpha / (h[1] * h[1]) // diffusion
          ;
      if (j > 0u)
        m_.add(idF, idF - mesh_.n_[0], valueBottomF);
      else
      {
        m_.add(idF, idF + mesh_.n_[0], valueBottomF);
        bcs_[0].bottom().ghostValues.set(i, valueBottomF);
      }

      // right
      double const valueRightF = -alpha / (h[0] * h[0]) // diffusion
          ;
      if (i < mesh_.n_[0] - 1)
        m_.add(idF, idF + 1, valueRightF);
      else
      {
        m_.add(idF, idF - 1, valueRightF);
        bcs_[0].right().ghostValues.set(j, valueRightF);
      }

      // top
      double const valueTopF = -alpha / (h[1] * h[1]) // diffusion
          ;
      if (j < mesh_.n_[1] - 1)
        m_.add(idF, idF + mesh_.n_[0], valueTopF);
      else
      {
        m_.add(idF, idF - mesh_.n_[0], valueTopF);
        bcs_[0].top().ghostValues.set(i, valueTopF);
      }

      // left
      double const valueLeftF = -alpha / (h[0] * h[0]) // diffusion
          ;
      if (i > 0u)
        m_.add(idF, idF - 1, valueLeftF);
      else
      {
        m_.add(idF, idF + 1, valueLeftF);
        bcs_[0].left().ghostValues.set(j, valueLeftF);
      }

      // coupling with adjoint
      if (pt[0] > control[0] - eps && pt[0] < control[1] + eps)
        m_.add(idF, idA, 1.0 / beta);

      // rhs
      rhs_.add(
          idF,
          uOld_[idF] / dt_ // time derivative
      );

      // adjoint
      // middle
      double const valueA =
          1. / dt_                                                 // time derivative
          + 2. * alpha * (1. / (h[0] * h[0]) + 1. / (h[1] * h[1])) // diffusion
          ;
      m_.add(idA, idA, valueA);

      // bottom
      double const valueBottomA = -alpha / (h[1] * h[1]) // diffusion
          ;
      if (j > 0u)
        m_.add(idA, idA - mesh_.n_[0], valueBottomA);
      else
      {
        m_.add(idA, idA + mesh_.n_[0], valueBottomA);
        bcs_[1].bottom().ghostValues.set(i, valueBottomA);
      }

      // right
      double const valueRightA = -alpha / (h[0] * h[0]) // diffusion
          ;
      if (i < mesh_.n_[0] - 1)
        m_.add(idA, idA + 1, valueRightA);
      else
      {
        m_.add(idA, idA - 1, valueRightA);
        bcs_[1].right().ghostValues.set(j, valueRightA);
      }

      // top
      double const valueTopA = -alpha / (h[1] * h[1]) // diffusion
          ;
      if (j < mesh_.n_[1] - 1)
        m_.add(idA, idA + mesh_.n_[0], valueTopA);
      else
      {
        m_.add(idA, idA - mesh_.n_[0], valueTopA);
        bcs_[1].top().ghostValues.set(i, valueTopA);
      }

      // left
      double const valueLeftA = -alpha / (h[0] * h[0]) // diffusion
          ;
      if (i > 0u)
        m_.add(idA, idA - 1, valueLeftA);
      else
      {
        m_.add(idA, idA + 1, valueLeftA);
        bcs_[1].left().ghostValues.set(j, valueLeftA);
      }

      // coupling with forward problem
      if (pt[0] > target[1] - eps && pt[0] < target[2] + eps)
      {
        m_.add(idA, idF, -1.0);
        rhs_.add(idA, -target[0]);
      }

      // rhs
      rhs_.add(
          idA,
          uOld_[idA] / dt_ // time derivative
      );
    }
  m_.close();
}

void ProblemFD2D::assemblyHeatCoupled()
{
  // // std::vector<double> uExt(n_, 2.0);
  // std::vector<double> uExt(n_);
  // double const * dataPtr = getField(nameExt_)->dataPtr();
  // std::copy(dataPtr, dataPtr + n_, uExt.data());

  // double const kAmpli = 10.;
  // for (uint k = 1u; k < n_ - 1; k++)
  // {
  //   // matrix
  //   m_.diag[k] = 1. / dt_                  // time
  //                + alpha_ * 2. / (h_ * h_) // diffusion
  //                + kAmpli                  // feedback control
  //       ;
  //   m_.diagUp[k] = -alpha_ / (h_ * h_);   // diffusion
  //   m_.diagDown[k] = -alpha_ / (h_ * h_); // diffusion

  //   // rhs
  //   rhs_[k] = uOld_[k] / dt_     // time
  //             + q_[k]            // source
  //             + kAmpli * uExt[k] // feedback control
  //       ;
  // }
  // m_.close();
}

void ProblemFD2D::print()
{
  if (it % printStep_ == 0)
  {
    for (uint v = 0u; v < nVars_; v++)
    {
      auto const filename =
          fmt::format("{}.{}.csv", (outputPrefix_ / varNames_[v]).string(), it);
      std::FILE * out = std::fopen(filename.c_str(), "w");
      fmt::println(
          out, "# problemfd2d csv, mesh size: {} {}", mesh_.n_[0] - 1, mesh_.n_[1] - 1);
      fmt::print(out, "{}, ", varNames_[v]);
      for (auto const & [fieldName, _]: fields_)
        fmt::print(out, "{}, ", fieldName);
      fmt::println(out, "x, y, z");
      for (uint j = 0; j < mesh_.n_[1]; j++)
        for (uint i = 0; i < mesh_.n_[0]; i++)
        {
          auto const id = i + j * mesh_.n_[0];
          auto const pt = mesh_.pt({i, j});
          fmt::print(out, "{:.6e},", u_[id + v * mesh_.nPts()]);
          for (auto const & [_, field]: fields_)
            fmt::print(out, "{:.6e}, ", field[id]);
          fmt::println(out, "{:.6e}, {:.6e}, 0.0", pt[0], pt[1]);
        }
      std::fclose(out);
    }
  }
}

void ProblemFD2D::printFields()
{
  auto const mesh = initMeshCouplingVolume(COUPLING_TYPE::MEDCOUPLING);
  for (auto const & [name, field]: fields_)
  {
    auto fieldMED = FieldCoupling::build(COUPLING_TYPE::MEDCOUPLING);
    fieldMED->init(
        name, mesh.get(), SUPPORT_TYPE::ON_NODES, NATURE_TYPE::INTENSIVE_MAXIMUM);
    fieldMED->setValues(field.data_);
    fieldMED->initIO(outputPrefix_);
    fieldMED->printVTK(0.0, 0);
  }
}

void ProblemFD2D::printSetup(std::string_view filename)
{
  std::abort();
  std::FILE * out = std::fopen(filename.data(), "w");

  fmt::println(out, "name: {}", name_);
  fmt::println(out, "debug: {}", debug_);
  fmt::println(out, "start: {:.6e} {:.6e}", mesh_.start_[0], mesh_.start_[1]);
  fmt::println(out, "end: {:.6e} {:.6e}", mesh_.end()[0], mesh_.end()[1]);
  fmt::println(out, "n_elems: {} {}", mesh_.n_[0] - 1, mesh_.n_[1] - 1);
  fmt::println(out, "n_vars: {}", nVars_);
  fmt::print(out, "var_names: ");
  for (auto const & varName: varNames_)
    fmt::print(out, "{} ", varName);
  fmt::print(out, "initial_value: ");
  for (uint v = 0u; v < nVars_; v++)
    fmt::print("0.0 ");
  fmt::println(out, "");
  for (auto const & [name, param]: params_.data_)
    fmt::println(out, "params: {}", name);
  fmt::println(out, "");
  fmt::println(out, "");
  fmt::println(out, "");

  std::fclose(out);
}

Marker ProblemFD2D::findRegion(std::string_view name)
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

std::unordered_map<FD_SOLVER_TYPE, Solver_T<ProblemFD2D::Matrix_T>>
    ProblemFD2D::solvers_ = {
        {FD_SOLVER_TYPE::JACOBI, &solveJacobi<ProblemFD2D::Matrix_T>},
        {FD_SOLVER_TYPE::GAUSS_SEIDEL, &solveGaussSeidel<ProblemFD2D::Matrix_T>},
        {FD_SOLVER_TYPE::CG, &solveConjugateGradient<ProblemFD2D::Matrix_T>},
        {FD_SOLVER_TYPE::BICGSTAB, &solveBiCGStab<ProblemFD2D::Matrix_T>},
        {FD_SOLVER_TYPE::VANKA2DCB, &solveVanka2DCB<ProblemFD2D::Matrix_T>},
        {FD_SOLVER_TYPE::VANKA2DSCI, &solveVanka2DSCI<ProblemFD2D::Matrix_T>},
};

} // namespace cocoa
