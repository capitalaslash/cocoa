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

ProblemFD2D::ProblemFD2D(): Problem{PROBLEM_TYPE::FD2D, COUPLING_TYPE::NONE}
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
  std::filesystem::path configFile = configs.at("config_file");
  std::ifstream in(configFile, std::ios::in);
  if (!in)
  {
    fmt::print(stderr, "configuration file {} not found!\n", configFile.string());
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
      // coupling
      else if (token == "coupling_type:")
      {
        bufferStream >> token;
        couplingType_ = str2coupling(token);
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
            fmt::print(stderr, "param type for {} not recognized\n", name);
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
        fmt::print(stderr, "key {} invalid\n", token);
        bufferStream >> token;
      }
    }
  }
  fmt::print("{} - equation type: {}\n", name_, eqn2str(eqnType_));
  assert(eqnType_ == EQN_TYPE::NONE || assemblies_.contains(eqnType_));

  fmt::print("parameters: {}\n", params_);

  // mesh
  mesh_.init(start, end, nElems);
  initMeshCoupling();

  // fields
  u_.resize(mesh_.nPts() * nVars_);
  uOld_.resize(mesh_.nPts() * nVars_);
  q_.resize(mesh_.nPts() * nVars_);
  c_[0].resize(mesh_.nPts() * nVars_);
  c_[1].resize(mesh_.nPts() * nVars_);
  std::vector<double> cValues(nVars_ * 2u, 0.0);
  if (params_.data_.contains("c"))
  {
    cValues = params_.get<FD_PARAM_TYPE::VECTOR>("c");
    assert(cValues.size() == nVars_ * 2u);
  }
  for (uint v = 0u; v < nVars_; v++)
  {
    u_.setRange(0u + v * mesh_.nPts(), mesh_.nPts() + v * mesh_.nPts(), uInit[v]);
    uOld_.setRange(0u + v * mesh_.nPts(), mesh_.nPts() + v * mesh_.nPts(), uInit[v]);
    q_.setRange(0u + v * mesh_.nPts(), mesh_.nPts() + v * mesh_.nPts(), qValue[v]);
    c_[0].setRange(
        0u + v * mesh_.nPts(), mesh_.nPts() + v * mesh_.nPts(), cValues[2 * v + 0]);
    c_[1].setRange(
        0u + v * mesh_.nPts(), mesh_.nPts() + v * mesh_.nPts(), cValues[2 * v + 1]);

    for (uint s = 0u; s < 4u; s++)
      bcs_[v].data_[s] =
          FDBC(side2D(s), bcTypes[v][s], bcValues[v][s], mesh_.n_[1 - s / 2]);
  }
  initFieldCoupling();

  // linear algebra
  m_.init(mesh_.nPts() * nVars_, nnz);
  rhs_.resize(mesh_.nPts() * nVars_);

  // io
  if (cleanOutput_)
    for (const auto & entry: std::filesystem::directory_iterator(outputPrefix_))
      std::filesystem::remove_all(entry.path());
  std::filesystem::create_directories(outputPrefix_);
}

void ProblemFD2D::initMeshCoupling()
{
  // coords format: x_0, y_0, z_0, x_1, ...
  std::vector<double> coords(mesh_.nPts() * 3);
  for (uint j = 0; j < mesh_.n_[1]; j++)
    for (uint i = 0; i < mesh_.n_[0]; i++)
    {
      uint const id = j * mesh_.n_[0] + i;
      auto const pt = mesh_.pt({i, j});
      coords[3 * id + 0U] = pt[0];
      coords[3 * id + 1U] = pt[1];
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

  meshCoupling_ = MeshCoupling::build(couplingType_);
  meshCoupling_->init("mesh_fd2d", 2U, 2U, coords, conn, offsets);
  // meshCoupling_->printVTK(outputPrefix_ / "mesh_fd2d");
}

void ProblemFD2D::initFieldCoupling()
{
  for (uint v = 0u; v < nVars_; v++)
  {
    auto [kvPairU, successU] =
        fieldsCoupling_.emplace(varNames_[v], FieldCoupling::build(couplingType_));
    assert(successU);
    kvPairU->second->init(varNames_[v], meshCoupling_.get(), SUPPORT_TYPE::ON_NODES);
    kvPairU->second->setValues(
        {u_.data() + v * mesh_.nPts(), u_.data() + (v + 1) * mesh_.nPts()}, 1u);
    kvPairU->second->initIO(outputPrefix_);

    auto const nameExt = (nVars_ > 1u) ? fmt::format("{}{}", nameExt_, v) : nameExt_;
    auto [kvPairExt, successExt] =
        fieldsCoupling_.emplace(nameExt, FieldCoupling::build(couplingType_));
    assert(successExt);
    kvPairExt->second->init(nameExt, meshCoupling_.get(), SUPPORT_TYPE::ON_NODES);
    kvPairExt->second->setValues(0.0, mesh_.nPts(), 1u);
  }
}

bool ProblemFD2D::run() { return time < finalTime_; }

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
// TODO: enum for sides
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

// static constexpr int sideOffset(std::array<uint, 2U> const & n, uint k)
// {
//   switch (k)
//   {
//   case 0U: // bottom
//     return n[0];
//   case 1U: // right
//     return -1;
//   case 2U: // top
//     return -n[0];
//   case 3U: // left
//     return 1;
//   default:
//     std::abort();
//   }
//   return 0;
// }

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
//   case 1U: // bottom-right
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

// static constexpr std::pair<uint, uint> cornerEnd(std::array<uint, 2U> const & n, uint
// k)
// {
//   switch (k)
//   {
//   case 0U:
//     return {0, 0}; // bottom left
//   case 1U:
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
  fmt::print("\n===\n");
  fmt::print("{}, time = {:.6e}, dt = {:.6e}\n", name_, time, dt_);

  // TODO: improve CFL evaluation by using better estimation of cell diameter
  if (computeCFL_)
  {
    double maxCFL = 0.0;
    for (uint k = 0U; k < u_.size(); k++)
    {
      double const cLocal = std::sqrt(c_[0][k] * c_[0][k] + c_[1][k] * c_[1][k]);
      maxCFL = std::max(cLocal * dt_ / std::min(mesh_.h_[0], mesh_.h_[1]), maxCFL);
    }
    fmt::print("maxCFL: {:.6e}\n", maxCFL);
  }
  // update
  uOld_ = u_;

  // assembly
  assemblies_.at(eqnType_)(this);

  // fmt::print("bc b: {}\n", bcs_[0].bottom());
  // fmt::print("bc r: {}\n", bcs_[0].right());
  // fmt::print("bc t: {}\n", bcs_[0].top());
  // fmt::print("bc l: {}\n", bcs_[0].left());

  std::array<double, 4U> const hSide = {
      mesh_.h_[1], mesh_.h_[1], mesh_.h_[0], mesh_.h_[0]};
  for (uint v = 0u; v < nVars_; v++)
  {
    // bc: Neumann sides
    for (uint s = 0u; s < 4u; s++)
    {
      auto const & bc = bcs_[v].data_[s];
      auto const dofList = sideDOF(mesh_.n_, side2D(s));
      // auto const offset = sideOffset(n_, k);

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
  if (debug_)
    m_.print_sparsity_pattern("fd2d_mat.dat");

  // solve
  auto const [numIters, residual] =
      solvers_.at(solverType_)(m_, rhs_, u_, tol_, maxIters_);
  fmt::print("num iters: {:4d}, ", numIters);
  fmt::print("relative residual: {:.8e}\n", residual);

  if (debug_)
  {
    fmt::print("matrix: {}\n", m_);
    fmt::print("rhs: {}\n", rhs_);
    fmt::print("sol: {}\n", u_);
  }
  // clean up
  m_.clear();
  rhs_.zero();

  // update coupling field
  for (uint v = 0u; v < nVars_; v++)
    getField(varNames_[v])
        ->setValues({u_.data() + v * mesh_.nPts(), u_.data() + (v + 1) * mesh_.nPts()});

  return numIters;
}

void ProblemFD2D::assemblyHeat()
{
  auto const h = mesh_.h_;
  auto const alpha = params_.get<FD_PARAM_TYPE::SCALAR>("alpha");

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
      double const valueBottom = -alpha / (h[1] * h[1])   // diffusion
                                 - 0.5 * c_[1][id] / h[1] // advection
          ;
      if (j > 0u)
        m_.add(id, id - mesh_.n_[0], valueBottom);
      else
      {
        m_.add(id, id + mesh_.n_[0], valueBottom);
        bcs_[0].bottom().ghostValues.set(i, valueBottom);
      }

      // right
      double const valueRight = -alpha / (h[0] * h[0])   // diffusion
                                + 0.5 * c_[0][id] / h[0] // advection
          ;
      if (i < mesh_.n_[0] - 1)
        m_.add(id, id + 1, valueRight);
      else
      {
        m_.add(id, id - 1, valueRight);
        bcs_[0].right().ghostValues.set(j, valueRight);
      }

      // top
      double const valueTop = -alpha / (h[1] * h[1])   // diffusion
                              + 0.5 * c_[1][id] / h[1] // advection
          ;
      if (j < mesh_.n_[1] - 1)
        m_.add(id, id + mesh_.n_[0], valueTop);
      else
      {
        m_.add(id, id - mesh_.n_[0], valueTop);
        bcs_[0].top().ghostValues.set(i, valueTop);
      }

      // left
      double const valueLeft = -alpha / (h[0] * h[0])   // diffusion
                               - 0.5 * c_[0][id] / h[0] // advection
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
  // for (uint k = 1U; k < n_ - 1; k++)
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
          fmt::format("{}.{}.dat", (outputPrefix_ / varNames_[v]).string(), it);
      std::FILE * out = std::fopen(filename.c_str(), "w");
      for (uint i = 0; i < mesh_.n_[0]; i++)
        for (uint j = 0; j < mesh_.n_[1]; j++)
        {
          auto const id = i + j * mesh_.n_[0];
          auto const pt = mesh_.pt({i, j});
          fmt::print(out, "{:.6e} {:.6e} {:.6e}\n", pt[0], pt[1], u_[id]);
        }
      std::fclose(out);

      getField(varNames_[v])->printVTK(time, it);
    }
  }
}

void ProblemFD2D::printSetup(std::string_view filename)
{
  std::abort();
  std::FILE * out = std::fopen(filename.data(), "w");

  fmt::print(out, "name: {}\n", name_);
  fmt::print(out, "debug: {}\n", debug_);
  fmt::print(out, "start: {:.6e} {:.6e}\n", mesh_.start_[0], mesh_.start_[1]);
  fmt::print(out, "end: {:.6e} {:.6e}\n", mesh_.end()[0], mesh_.end()[1]);
  fmt::print(out, "n_elems: {} {}\n", mesh_.n_[0] - 1, mesh_.n_[1] - 1);
  fmt::print(out, "coupling_type: {}\n", couplingType2str(couplingType_));
  fmt::print(out, "n_vars: {}\nvar_names: ", nVars_);
  for (auto const & varName: varNames_)
    fmt::print(out, "{} ", varName);
  fmt::print(out, "initial_value: ");
  for (uint v = 0u; v < nVars_; v++)
    fmt::print("0.0 ");
  fmt::print(out, "\n");
  for (auto const & [name, param]: params_.data_)
    fmt::print(out, "params: {}\n", name);
  fmt::print(out, "\n");
  fmt::print(out, "\n");
  fmt::print(out, "\n");

  std::fclose(out);
}

std::unordered_map<FD_SOLVER_TYPE, Solver_T<ProblemFD2D::Matrix_T>>
    ProblemFD2D::solvers_ = {
        {FD_SOLVER_TYPE::GAUSS_SEIDEL, &solveGaussSeidel<ProblemFD2D::Matrix_T>},
        {FD_SOLVER_TYPE::JACOBI, &solveJacobi<ProblemFD2D::Matrix_T>},
        {FD_SOLVER_TYPE::VANKA2DCB, &solveVanka2DCB<ProblemFD2D::Matrix_T>},
        {FD_SOLVER_TYPE::VANKA2DSCI, &solveVanka2DSCI<ProblemFD2D::Matrix_T>},
};
