#include "problem/problem_fd2d.hpp"

// std
#include <cassert>
#include <cmath>
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

void ProblemFD2D::setup(Problem::ParamList_T const & params)
{
  // default values
  name_ = "empty";
  // mesh
  Vec2D_T start = {0.0, 0.0};
  Vec2D_T end = {1.0, 1.0};
  std::array<uint, 2U> nElems = {10U, 10U};
  // fields
  double uInit = 0.0;
  double qValue = 1.0;
  // physical constants
  alpha_ = 0.2;
  std::array<double, 2U> cValues{0.0, 0.0};
  // bcs
  std::array<double, 4U> bcValues{0.0, 0.0, 0.0, 0.0};
  // time
  time = 0.0;
  finalTime_ = 1.0;
  dt_ = 0.1;
  // linear algebra
  maxIters_ = 1000U;
  toll_ = 1.e-6;

  // read configuration from file
  std::filesystem::path configFile = params.at("config_file");
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
      else if (token == "coupling_type:")
      {
        bufferStream >> token;
        couplingType_ = str2coupling(token);
      }
      else if (token == "var_name:")
        bufferStream >> varName_;
      else if (token == "initial_value:")
        bufferStream >> uInit;
      else if (token == "q:")
        bufferStream >> qValue;
      else if (token == "alpha:")
        bufferStream >> alpha_;
      else if (token == "c:")
      {
        bufferStream >> cValues[0];
        bufferStream >> cValues[1];
      }
      else if (token == "start_time:")
        bufferStream >> time;
      else if (token == "final_time:")
        bufferStream >> finalTime_;
      else if (token == "dt:")
        bufferStream >> dt_;
      else if (token == "assembly_name:")
      {
        bufferStream >> token;
        eqnType_ = str2eqn(token);
      }
      else if (token == "solver_type:")
      {
        bufferStream >> token;
        solverType_ = str2fdsolver(token);
      }
      else if (token == "max_iters:")
        bufferStream >> maxIters_;
      else if (token == "toll:")
        bufferStream >> toll_;
      else if (token == "bcs:")
      {
        for (uint k = 0U; k < 4U; k++)
        {
          bufferStream >> token;
          bcs_[k].type = str2fdbc(token);
          bufferStream >> bcValues[k];
        }
      }
      else if (token == "output_prefix:")
        bufferStream >> outputPrefix_;
      else
      {
        fmt::print(stderr, "key {} invalid\n", token);
        bufferStream >> token;
      }
    }
  }
  fmt::print("{} - equation type: {}\n", name_, eqn2str(eqnType_));
  assert(eqnType_ == EQN_TYPE::NONE || assemblies_.contains(eqnType_));

  // mesh
  start_ = start;
  h_ = {(end[0] - start[0]) / nElems[0], (end[1] - start[1]) / nElems[1]};
  n_ = {nElems[0] + 1, nElems[1] + 1};
  initMeshCoupling();

  // fields
  u_.resize(n_[0] * n_[1], uInit);
  uOld_.resize(n_[0] * n_[1], uInit);
  q_.resize(n_[0] * n_[1], qValue);
  initFieldCoupling();

  // physical properties
  c_[0].resize(n_[0] * n_[1], cValues[0]);
  c_[1].resize(n_[0] * n_[1], cValues[1]);

  // bcs
  for (uint k = 0U; k < 4U; k++)
    bcs_[k].values.resize(n_[k % 2], bcValues[k]);

  // linear algebra
  m_.init(n_[0] * n_[1]);
  rhs_.resize(n_[0] * n_[1]);

  // io
  std::filesystem::create_directories(outputPrefix_);
}

void ProblemFD2D::initMeshCoupling()
{
  // coords format: x_0, y_0, z_0, x_1, ...
  std::vector<double> coords(n_[0] * n_[1] * 3);
  for (uint j = 0; j < n_[1]; j++)
    for (uint i = 0; i < n_[0]; i++)
    {
      uint const id = j * n_[0] + i;
      coords[3 * id + 0U] = start_[0] + i * h_[0];
      coords[3 * id + 1U] = start_[1] + j * h_[1];
      coords[3 * id + 2U] = 0.0;
    }

  // conn format: elem0_numpts, id_0, id_1, ..., elem1_numpts, ...
  auto const nElems = (n_[0] - 1) * (n_[1] - 1);
  std::vector<uint> conn(nElems * (1 + 4));
  auto elemCount = 0U;
  for (uint j = 0; j < n_[1] - 1; j++)
    for (uint i = 0; i < n_[0] - 1; i++)
    {
      uint const id = j * n_[0] + i;
      conn[(1 + 4) * elemCount] = 4; // MEDCellTypeToIKCell(MED_CELL_TYPE::QUAD4);
      conn[(1 + 4) * elemCount + 1] = id;
      conn[(1 + 4) * elemCount + 2] = id + 1;
      conn[(1 + 4) * elemCount + 3] = id + n_[0] + 1;
      conn[(1 + 4) * elemCount + 4] = id + n_[0];
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
  auto [kvPairU, successU] =
      fieldsCoupling_.emplace(varName_, FieldCoupling::build(couplingType_));
  assert(successU);
  kvPairU->second->init(varName_, meshCoupling_.get(), SUPPORT_TYPE::ON_NODES);
  kvPairU->second->setValues(u_);
  kvPairU->second->initIO(outputPrefix_);

  auto [kvPairExt, successExt] =
      fieldsCoupling_.emplace(nameExt_, FieldCoupling::build(couplingType_));
  assert(successExt);
  kvPairExt->second->init(nameExt_, meshCoupling_.get(), SUPPORT_TYPE::ON_NODES);
  kvPairExt->second->setValues(0.0, u_.size());
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
const std::vector<uint> sideDOF(std::array<uint, 2U> const & n, uint const side)
{
  uint const dofSize = n[side % 2];
  std::vector<uint> dofList(dofSize);

  switch (side)
  {
  case 0U: // bottom
  {
    std::iota(dofList.begin(), dofList.end(), 0U);
    break;
  }
  case 1U: // right
  {
    for (uint k = 0U; k < dofSize; k++)
      dofList[k] = (k + 1) * n[0] - 1;
    break;
  }
  case 2U: // top
  {
    std::iota(dofList.begin(), dofList.end(), n[0] * (n[1] - 1));
    break;
  }
  case 3U: // left
  {
    for (uint k = 0U; k < dofSize; k++)
      dofList[k] = k * n[0];
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

static constexpr std::array<uint, 4U> cornerDOF(std::array<uint, 2U> const & n)
{
  return std::array<uint, 4U>{{
      0U,                // bottom-left
      n[0] - 1,          // bottom-right
      n[0] * n[1] - 1,   // top-right
      n[0] * (n[1] - 1), // top-left
  }};
}

static constexpr std::array<std::array<uint, 2U>, 4U> cornerSides = {{
    {{0, 3}}, // bottom-left
    {{1, 0}}, // bottom-right
    {{2, 1}}, // top-right
    {{3, 2}}, // top-left
}};

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

static constexpr std::pair<uint, uint> cornerEnd(std::array<uint, 2U> const & n, uint k)
{
  switch (k)
  {
  case 0U:
    return {0, 0}; // bottom left
  case 1U:
    return {0, n[0] - 1}; // bottom right
  case 2U:
    return {n[0] - 1, n[1] - 1}; // top left
  case 3U:
    return {n[1] - 1, 0}; // top right
  default:
    std::abort();
  }
  return {0, 0};
}

// TODO: unify with 1D version by using templates
double computeResidual(
    MatrixCSR const & m,
    std::vector<double> const & x,
    std::vector<double> const & b,
    double const area = 1.0)
{
  std::vector<double> res(x.size());
  auto const tmp = m * x;
  for (uint k = 0; k < res.size(); k++)
    res[k] = (b[k] - tmp[k]);
  double const resNorm = std::sqrt(norm2sq(res) * area);
  return resNorm;
}

void ProblemFD2D::solve()
{
  fmt::print("\n===\n");
  fmt::print("{}, time = {:.6e}, dt = {:.6e}\n", name_, time, dt_);

  // TODO: improve CFL evaluation by using better estimation of cell diameter
  double maxCFL = 0.0;
  for (uint k = 0U; k < u_.size(); k++)
  {
    double const cLocal = std::sqrt(c_[0][k] * c_[0][k] + c_[1][k] * c_[1][k]);
    maxCFL = std::max(cLocal * dt_ / std::min(h_[0], h_[1]), maxCFL);
  }
  fmt::print("maxCFL: {:.6e}\n", maxCFL);

  // update
  for (uint k = 0; k < u_.size(); k++)
    uOld_[k] = u_[k];

  // assembly
  assemblies_.at(eqnType_)(this);

  // TODO: decide if sign should mean always entrant/always outgoing (+1/-1 always) or
  // follow direction axis (current)
  std::array<double, 4U> const sideSign = {-1.0, 1.0, 1.0, -1.0};
  std::array<double, 4U> const hSide = {h_[1], h_[0], h_[1], h_[0]};

  // Neumann sides
  for (uint s = 0U; s < 4U; s++)
  {
    auto const & bc = bcs_[s];
    auto const dofList = sideDOF(n_, s);
    // auto const offset = sideOffset(n_, k);

    if (bc.type == FD_BC_TYPE::NEUMANN)
    {
      for (uint k = 0U; k < dofList.size(); k++)
      {
        uint const dof = dofList[k];
        rhs_[dof] += sideSign[s] * 2.0 * alpha_ * bc.values[k] / hSide[s];
      }
    }
  }

  // Dirichlet sides
  for (uint s = 0U; s < 4U; s++)
  {
    auto const & bc = bcs_[s];

    if (bc.type == FD_BC_TYPE::DIRICHLET)
    {
      auto const dofList = sideDOF(n_, s);
      for (uint k = 0U; k < dofList.size(); k++)
      {
        uint const dof = dofList[k];
        m_.clearRow(dof);
        m_.add(dof, dof, 1.0);
        rhs_[dof] = bc.values[k];
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
  //         bcs_[cornerSides[k][0]].values[] * h_[1] + bcs_[cornerSides[k][1]].values[]
  //         * h_[0];
  //   }
  // }

  m_.close();

  // solve
  auto const [numIters, residual] = solvers_.at(solverType_)(this);
  fmt::print("num iters: {:4d}, ", numIters);
  fmt::print("relative residual: {:.8e}\n", residual);

  // fmt::print("matrix: {}\n", m_);
  // fmt::print("rhs: {}\n", rhs_);
  // fmt::print("sol: {}\n", u_);

  // clean up
  m_.clear();
  for (uint k = 0; k < rhs_.size(); k++)
    rhs_[k] = 0.0;

  // update coupling field
  getField(varName_)->setValues(u_);
}

void ProblemFD2D::assemblyHeat()
{
  for (uint j = 0U; j < n_[1]; j++)
    for (uint i = 0U; i < n_[0]; i++)
    {
      auto const id = i + j * n_[0];

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
          1. / dt_ // time derivative
          + 2. * alpha_ * (1. / (h_[0] * h_[0]) + 1. / (h_[1] * h_[1])) // diffusion
          ;
      m_.add(id, id, value);

      // bottom
      uint const idBottom = (j != 0U) ? id - n_[0] : id + n_[0];
      double const valueBottom = -alpha_ / (h_[1] * h_[1]) // diffusion
                                 - 0.5 * c_[1][id] / h_[1] // advection
          ;
      m_.add(id, idBottom, valueBottom);

      // right
      uint const idRight = (i != n_[0] - 1) ? id + 1 : id - 1;
      double const valueRight = -alpha_ / (h_[0] * h_[0]) // diffusion
                                + 0.5 * c_[0][id] / h_[0] // advection
          ;
      m_.add(id, idRight, valueRight);

      // top
      uint const idTop = (j != n_[1] - 1) ? id + n_[0] : id - n_[0];
      double const valueTop = -alpha_ / (h_[1] * h_[1]) // diffusion
                              + 0.5 * c_[1][id] / h_[1] // advection
          ;
      m_.add(id, idTop, valueTop);

      // left
      uint const idLeft = (i != 0U) ? id - 1 : id + 1;
      double const valueLeft = -alpha_ / (h_[0] * h_[0]) // diffusion
                               - 0.5 * c_[0][id] / h_[0] // advection
          ;
      m_.add(id, idLeft, valueLeft);

      rhs_[id] = uOld_[id] / dt_ // time derivative
                 + q_[id]        // source
          ;
    }
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
}

std::pair<uint, double> ProblemFD2D::solveJacobi()
{
  std::vector<double> uNew(u_.size());
  double const rhsNorm = std::sqrt(norm2sq(rhs_) * h_[0] * h_[1]);
  fmt::print("rhsNorm: {:.8e}\n", rhsNorm);

  for (uint n = 0; n < maxIters_; n++)
  {
    double const resNorm = computeResidual(m_, u_, rhs_, h_[0] * h_[1]);
    // fmt::print("iter: {:3d}, current residual: {:.8e}\n", j, resNorm);
    if (resNorm < toll_ * rhsNorm)
    {
      return std::pair{n, resNorm / rhsNorm};
    }

    for (uint k = 0U; k < u_.size(); k++)
    {
      double valueNew = rhs_[k];
      double diag = 0.0;
      for (auto const [clm, value]: m_.data_[k])
      {
        if (clm == k)
        {
          diag = value;
        }
        else
        {
          valueNew -= value * u_[clm];
        }
      }
      uNew[k] = valueNew / diag;
    }

    // update current solution
    for (uint k = 0U; k < uNew.size(); k++)
      u_[k] = uNew[k];
  }

  return std::pair{maxIters_, computeResidual(m_, u_, rhs_, h_[0] * h_[1]) / rhsNorm};
}

std::pair<uint, double> ProblemFD2D::solveGaussSeidel()
{
  double const rhsNorm = std::sqrt(norm2sq(rhs_) * h_[0] * h_[1]);
  fmt::print("rhsNorm: {:.8e}\n", rhsNorm);

  for (uint n = 0; n < maxIters_; n++)
  {
    double const resNorm = computeResidual(m_, u_, rhs_, h_[0] * h_[1]);
    // fmt::print("iter: {:3d}, current residual: {:.8e}\n", j, resNorm);
    if (resNorm < toll_ * rhsNorm)
    {
      return std::pair{n, resNorm / rhsNorm};
    }

    for (uint k = 0U; k < u_.size(); k++)
    {
      double valueNew = rhs_[k];
      double diag = 0.0;
      for (auto const [clm, value]: m_.data_[k])
      {
        if (clm == k)
        {
          diag = value;
        }
        else
        {
          valueNew -= value * u_[clm];
        }
      }
      u_[k] = valueNew / diag;
    }
  }

  return std::pair{maxIters_, computeResidual(m_, u_, rhs_, h_[0] * h_[1]) / rhsNorm};
}

inline double solveLine(
    MatrixCSR const & m,
    std::vector<double> const & rhs,
    std::vector<double> const & u,
    size_t const id)
{
  assert(m.data_[id][0].clm == id);
  double value = rhs[id];
  for (uint j = 1; j < m.data_[id].size(); j++)
  {
    value -= m.data_[id][j].value * u[m.data_[id][j].clm];
  }
  return value / m.data_[id][0].value;
}

std::pair<uint, double> ProblemFD2D::solveVankaSCI()
{
  double const rhsNorm = std::sqrt(norm2sq(rhs_) * h_[0] * h_[1]);
  fmt::print("rhsNorm: {:.8e}\n", rhsNorm);

  for (uint n = 0U; n < maxIters_; n++)
  {
    double const resNorm = computeResidual(m_, u_, rhs_, h_[0] * h_[1]);

    // fmt::print("iter: {:3d}, current residual: {:.8e}\n", j, resNorm);
    if (resNorm < toll_ * rhsNorm)
    {
      return std::pair{n, resNorm / rhsNorm};
    }

    // for (uint k = 0U; k < u_.size(); k++)
    //   uOld_[k] = u_[k];

    // sides + corners + inside
    // sides
    for (uint k = 0U; k < 4U; k++)
    {
      auto const dofList = sideDOF(n_, k);
      for (auto const & dof: dofList)
        u_[dof] = solveLine(m_, rhs_, u_, dof);
    }

    // corners
    for (auto const & id: cornerDOF(n_))
      u_[id] = solveLine(m_, rhs_, u_, id);

    // inside
    for (uint j = 1U; j < n_[1] - 1; j += 1U)
      for (uint i = 1U; i < n_[0] - 1; i += 1U)
      {
        auto const id = i + j * n_[0];
        u_[id] = solveLine(m_, rhs_, u_, id);
      }
  }

  return std::pair{maxIters_, computeResidual(m_, u_, rhs_, h_[0] * h_[1]) / rhsNorm};
}

std::pair<uint, double> ProblemFD2D::solveVankaCB()
{
  double const rhsNorm = std::sqrt(norm2sq(rhs_) * h_[0] * h_[1]);
  fmt::print("rhsNorm: {:.8e}\n", rhsNorm);

  for (uint n = 0U; n < maxIters_; n++)
  {
    double const resNorm = computeResidual(m_, u_, rhs_, h_[0] * h_[1]);

    // fmt::print("iter: {:3d}, current residual: {:.8e}\n", j, resNorm);
    if (resNorm < toll_ * rhsNorm)
    {
      return std::pair{n, resNorm / rhsNorm};
    }

    // checkerboard
    for (uint k = 0U; k < n_[0] * n_[1]; k += 2U)
      u_[k] = solveLine(m_, rhs_, u_, k);
    for (uint k = 1U; k < n_[0] * n_[1]; k += 2U)
      u_[k] = solveLine(m_, rhs_, u_, k);
  }

  return std::pair{maxIters_, computeResidual(m_, u_, rhs_, h_[0] * h_[1]) / rhsNorm};
}

void ProblemFD2D::print()
{
  auto const filename =
      fmt::format("{}.{}.dat", (outputPrefix_ / varName_).string(), it);
  std::FILE * out = std::fopen(filename.c_str(), "w");
  for (uint i = 0; i < n_[0]; i++)
    for (uint j = 0; j < n_[1]; j++)
    {
      auto const id = i + j * n_[0];
      fmt::print(
          out,
          "{:.6e} {:.6e} {:.6e}\n",
          start_[0] + i * h_[0],
          start_[1] + j * h_[1],
          u_[id]);
    }
  std::fclose(out);

  getField(varName_)->printVTK(time, it);
}

std::unordered_map<EQN_TYPE, ProblemFD2D::Assembly_T> ProblemFD2D::assemblies_ = {
    {EQN_TYPE::HEAT, [](ProblemFD2D * p) { p->assemblyHeat(); }},
    // {EQN_TYPE::HEAT_COUPLED, [](ProblemFD2D * p) { p->assemblyHeatCoupled(); }},
};

std::unordered_map<FD_SOLVER_TYPE, ProblemFD2D::Solver_T> ProblemFD2D::solvers_ = {
    {FD_SOLVER_TYPE::JACOBI2D, [](ProblemFD2D * p) { return p->solveJacobi(); }},
    {FD_SOLVER_TYPE::GAUSSSEIDEL2D,
     [](ProblemFD2D * p) { return p->solveGaussSeidel(); }},
    {FD_SOLVER_TYPE::VANKA2DCB, [](ProblemFD2D * p) { return p->solveVankaCB(); }},
    {FD_SOLVER_TYPE::VANKA2DSCI, [](ProblemFD2D * p) { return p->solveVankaSCI(); }},
};
