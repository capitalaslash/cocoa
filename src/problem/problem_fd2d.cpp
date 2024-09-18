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
  // time
  time = 0.0;
  finalTime_ = 1.0;
  dt_ = 0.1;

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
      if (token == "name:")
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
      else if (token == "u_init:")
        bufferStream >> uInit;
      else if (token == "q:")
        bufferStream >> qValue;
      else if (token == "alpha:")
        bufferStream >> alpha_;
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
      else if (token == "bcs:")
      {
        for (uint k = 0U; k < 4U; k++)
        {
          bufferStream >> token;
          bcs_[k].type = str2fdbc(token);
          bufferStream >> bcs_[k].value;
        }
      }
      else if (token == "out_file:")
        bufferStream >> outFile_;
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
  auto const dirPath = std::filesystem::path(outFile_).parent_path();
  initMeshCoupling(dirPath / "mesh");

  // fields
  u_.resize(n_[0] * n_[1], uInit);
  uOld_.resize(n_[0] * n_[1], uInit);
  q_.resize(n_[0] * n_[1], qValue);
  initFieldCoupling(dirPath / "u");

  // linear algebra
  m_.init(n_[0] * n_[1]);
  rhs_.resize(n_[0] * n_[1]);

  // io
  std::filesystem::create_directories(dirPath);
}

void ProblemFD2D::initMeshCoupling(std::filesystem::path const & fileName)
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
  meshCoupling_->init(fileName.string(), 2U, 2U, coords, conn, offsets);
  meshCoupling_->printVTK();
}

void ProblemFD2D::initFieldCoupling(std::filesystem::path const & fileName)
{
  auto [kvPairU, successU] =
      fieldsCoupling_.emplace("u", FieldCoupling::build(couplingType_));
  assert(successU);
  kvPairU->second->init(
      fileName.filename().string(), meshCoupling_.get(), SUPPORT_TYPE::ON_NODES);
  kvPairU->second->setValues(u_);
  kvPairU->second->initIO(fileName.string() + "_med.");

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

// static constexpr std::vector<uint> // requires gcc >= 12
static const std::vector<uint> sideDOF(std::array<uint, 2U> const & n, uint const side)
{
  switch (side)
  {
  case 0U: // bottom
  {
    std::vector<uint> dofList(n[0] - 2);
    std::iota(dofList.begin(), dofList.end(), 1U);
    return dofList;
  }
  case 1U: // right
  {
    std::vector<uint> dofList(n[1] - 2);
    for (uint k = 0U; k < n[1] - 2; k++)
      dofList[k] = (k + 2) * n[0] - 1;
    return dofList;
  }
  case 2U: // top
  {
    std::vector<uint> dofList(n[0] - 2);
    std::iota(dofList.begin(), dofList.end(), n[0] * (n[1] - 1) + 1);
    return dofList;
  }
  case 3U: // left
  {
    std::vector<uint> dofList(n[1] - 2);
    for (uint k = 0U; k < n[1] - 2; k++)
      dofList[k] = (k + 1) * n[0];
    return dofList;
  }
  default:
    std::abort();
  }

  return std::vector<uint>(0U);
}

static constexpr std::array<uint, 4U> cornerDOF(std::array<uint, 2U> const & n)
{
  return std::array<uint, 4U>{{0U, n[0] - 1, n[0] * (n[1] - 1), n[0] * n[1] - 1}};
}

static constexpr int sideOffset(std::array<uint, 2U> const & n, uint k)
{
  switch (k)
  {
  case 0U:
    return n[0];
  case 1U:
    return -1;
  case 2U:
    return -n[0];
  case 3U:
    return 1;
  default:
    std::abort();
  }
  return 0;
}

static constexpr std::array<std::array<uint, 2U>, 4U> cornerSides = {{
    {{0, 3}},
    {{0, 1}},
    {{2, 3}},
    {{2, 1}},
}};

static constexpr int cornerOffset(std::array<uint, 2U> const & n, uint k)
{
  switch (k)
  {
  case 0U:
    return n[0] + 1; // bottom left
  case 1U:
    return n[0] - 1; // bottom right
  case 2U:
    return -n[0] + 1; // top left
  case 3U:
    return -n[0] - 1; // top right
  default:
    std::abort();
  }
  return 0;
}

void ProblemFD2D::solve()
{
  fmt::print("{}, time = {:.6e}, dt = {:.6e}\n", name_, time, dt_);
  Vec2D_T const h = {1. / (n_[0] - 1), 1. / (n_[0] - 1)};

  // update
  for (uint k = 0; k < u_.size(); k++)
    uOld_[k] = u_[k];

  // sides
  for (uint k = 0U; k < 4U; k++)
  {
    auto const & bc = bcs_[k];
    auto const dofList = sideDOF(n_, k);
    auto const offset = sideOffset(n_, k);
    switch (bc.type)
    {
    case FD_BC_TYPE::DIRICHLET:
    {
      for (auto const & dof: dofList)
      {
        m_.triplets_.emplace_back(dof, dof, 1.0);
        rhs_[dof] = bc.value;
      }
      break;
    }
    case FD_BC_TYPE::NEUMANN:
    {
      for (auto const & dof: dofList)
      {
        m_.triplets_.emplace_back(dof, dof, 1.0);
        m_.triplets_.emplace_back(dof, dof + offset, -1.0);
        rhs_[dof] = bc.value * h[1];
      }
      break;
    }
    default:
    {
      fmt::print(stderr, "bc {} not specified!\n", k);
      std::abort();
    }
    }
  }

  // corners
  for (uint k = 0; k < 4U; k++)
  {
    auto const dof = cornerDOF(n_)[k];
    m_.triplets_.emplace_back(dof, dof, 1.0);
    if (bcs_[cornerSides[k][0]].type == FD_BC_TYPE::DIRICHLET)
    {
      rhs_[dof] = bcs_[cornerSides[k][0]].value;
    }
    else if (bcs_[cornerSides[k][1]].type == FD_BC_TYPE::DIRICHLET)
    {
      rhs_[dof] = bcs_[cornerSides[k][1]].value;
    }
    else
    {
      m_.triplets_.emplace_back(dof, dof + cornerOffset(n_, k), -1.0);
      rhs_[dof] =
          bcs_[cornerSides[k][0]].value * h[1] + bcs_[cornerSides[k][1]].value * h[0];
    }
  }

  // fmt::print(stderr, "{}\n", m_.triplets_);

  // assembly
  assemblies_.at(eqnType_)(this);

  // // solve
  solvers_.at(solverType_)(this);

  // // residual
  std::vector<double> res(u_.size());
  auto const tmp = m_ * u_;
  for (uint k = 0; k < res.size(); k++)
    res[k] = rhs_[k] - tmp[k];
  double resNorm = 0.0;
  for (uint k = 0; k < res.size(); k++)
    resNorm += res[k] * res[k];
  fmt::print("residual norm: {:.6e}\n", std::sqrt(resNorm));

  m_.clear();
  for (uint k = 0; k < res.size(); k++)
    rhs_[k] = 0.0;

  getField("u")->setValues(u_);
}

void ProblemFD2D::assemblyHeat()
{
  for (uint j = 1; j < n_[1] - 1; j++)
    for (uint i = 1; i < n_[0] - 1; i++)
    {
      auto const id = i + j * n_[0];
      m_.triplets_.emplace_back(
          id, id, 1. / dt_ + alpha_ * (2. / (h_[0] * h_[0]) + 2. / (h_[1] * h_[1])));
      m_.triplets_.emplace_back(id, id - n_[0], -alpha_ / (h_[1] * h_[1]));
      m_.triplets_.emplace_back(id, id + 1, -alpha_ / (h_[0] * h_[0]));
      m_.triplets_.emplace_back(id, id + n_[0], -alpha_ / (h_[1] * h_[1]));
      m_.triplets_.emplace_back(id, id - 1, -alpha_ / (h_[0] * h_[0]));
      rhs_[id] = uOld_[id] / dt_ + q_[id];
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

void ProblemFD2D::solveVankaSCI()
{
  for (uint n = 0U; n < 20U; n++)
  {
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
}

void ProblemFD2D::solveVankaCB()
{
  for (uint n = 0U; n < 20U; n++)
  {
    // checkerboard
    for (uint k = 0U; k < n_[0] * n_[1]; k += 2U)
      u_[k] = solveLine(m_, rhs_, u_, k);
    for (uint k = 1U; k < n_[0] * n_[1]; k += 2U)
      u_[k] = solveLine(m_, rhs_, u_, k);
  }
}

void ProblemFD2D::print()
{
  std::FILE * out =
      std::fopen((outFile_ + "." + std::to_string(it) + ".dat").c_str(), "w");
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

  getField("u")->printVTK(time, it);
}

std::unordered_map<EQN_TYPE, ProblemFD2D::Assembly_T> ProblemFD2D::assemblies_ = {
    {EQN_TYPE::HEAT, [](ProblemFD2D * p) { p->assemblyHeat(); }},
    // {EQN_TYPE::HEAT_COUPLED, [](ProblemFD2D * p) { p->assemblyHeatCoupled(); }},
};

std::unordered_map<FD_SOLVER_TYPE, ProblemFD2D::Solver_T> ProblemFD2D::solvers_ = {
    {FD_SOLVER_TYPE::VANKA2DCB, [](ProblemFD2D * p) { p->solveVankaCB(); }},
    {FD_SOLVER_TYPE::VANKA2DSCI, [](ProblemFD2D * p) { p->solveVankaSCI(); }},
};
