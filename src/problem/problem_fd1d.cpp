#include "problem/problem_fd1d.hpp"

// std
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <unordered_map>

// libfmt
#include <fmt/core.h>

// local
#include "coupling/field_coupling.hpp"
#include "coupling/mesh_coupling.hpp"
#include "enums.hpp"
#include "problem/fdutils.hpp"

void ProblemFD1D::setup(Problem::ParamList_T const & params)
{
  // default values
  name_ = "empty";
  // mesh
  double start = 0.0;
  double end = 1.0;
  uint nElems = 10U;
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
  std::filesystem::path const configFile = params.at("config_file");
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
        bufferStream >> start;
      else if (token == "end:")
        bufferStream >> end;
      else if (token == "n_elems:")
        bufferStream >> nElems;
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
      else if (token == "bc_start:")
      {
        bufferStream >> token;
        bcStart_.type = str2fdbc(token);
        double value;
        bufferStream >> value;
        bcStart_.values.assign(1, value);
      }
      else if (token == "bc_end:")
      {
        bufferStream >> token;
        bcEnd_.type = str2fdbc(token);
        double value;
        bufferStream >> value;
        bcEnd_.values.assign(1, value);
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
  h_ = (end - start) / nElems;
  n_ = nElems + 1;
  auto const dirPath = std::filesystem::path(outFile_).parent_path();
  initMeshCoupling(dirPath / "mesh");

  // fields
  u_.resize(n_, uInit);
  uOld_.resize(n_, uInit);
  q_.resize(n_, qValue);
  initFieldCoupling(dirPath / "u");

  // linear algebra
  m_.init(n_);
  rhs_.resize(n_);

  // io
  std::filesystem::create_directories(dirPath);
}

void ProblemFD1D::initMeshCoupling(std::filesystem::path const & fileName)
{
  // coords format: x_0, y_0, z_0, x_1, ...
  std::vector<double> coords(n_ * 3);
  for (uint k = 0; k < n_; k++)
  {
    coords[3 * k] = start_ + k * h_;
    coords[3 * k + 1] = 0.0;
    coords[3 * k + 2] = 0.0;
  }

  // conn format: elem0_numpts, id_0, id_1, ..., elem1_numpts, ...
  auto const nElems = n_ - 1;
  std::vector<uint> conn(nElems * 3);
  for (uint k = 0; k < nElems; k++)
  {
    conn[3 * k] = MEDCellTypeToIKCell(MED_CELL_TYPE::LINE2);
    conn[3 * k + 1] = k;
    conn[3 * k + 2] = k + 1;
  }

  // offsets format: sum_0^k elemk_numpts + 1,
  std::vector<uint> offsets(nElems + 1);
  offsets[0] = 0;
  for (uint k = 0; k < nElems; k++)
  {
    offsets[k + 1] = offsets[k] + 3;
  }

  meshCoupling_ = MeshCoupling::build(couplingType_);
  meshCoupling_->init(fileName.string(), 1U, 1U, coords, conn, offsets);
}

void ProblemFD1D::initFieldCoupling(std::filesystem::path const & fileName)
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

bool ProblemFD1D::run() { return time < finalTime_; }

void ProblemFD1D::advance()
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

void ProblemFD1D::solve()
{
  fmt::print("\n===\n");
  fmt::print("{}, time = {:.6e}, dt = {:.6e}\n", name_, time, dt_);

  // update
  for (uint k = 0; k < u_.size(); k++)
    uOld_[k] = u_[k];

  // assembly
  assemblies_.at(eqnType_)(this);

  // bc start
  switch (bcStart_.type)
  {
  case FD_BC_TYPE::DIRICHLET:
  {
    m_.clearRow(0U);
    m_.diag()[0] = 1.0;
    rhs_[0] = bcStart_.values[0];
    break;
  }
  case FD_BC_TYPE::NEUMANN:
  {
    rhs_[0] += -2.0 * alpha_ * bcStart_.values[0] / h_;
    break;
  }
  default:
  {
    fmt::print(stderr, "no bc start specified!\n");
    std::abort();
  }
  }

  // bc end
  switch (bcEnd_.type)
  {
  case FD_BC_TYPE::DIRICHLET:
  {
    m_.clearRow(n_ - 1);
    m_.diag()[n_ - 1] = 1.0;
    rhs_[n_ - 1] = bcEnd_.values[0];
    break;
  }
  case FD_BC_TYPE::NEUMANN:
  {
    rhs_[n_ - 1] += 2.0 * alpha_ * bcEnd_.values[0] / h_;
    break;
  }
  default:
  {
    fmt::print(stderr, "no bc end specified!\n");
    std::abort();
  }
  }

  double const rhsNorm = std::sqrt(norm2sq(rhs_) * h_);
  fmt::print("rhsNorm: {:.8e}\n", rhsNorm);

  // solve
  auto const [numIters, residual] = solvers_.at(solverType_)(this);
  fmt::print("num iters: {:4d}, ", numIters);
  fmt::print("relative residual: {:.8e}\n", residual / rhsNorm);

  // clean up
  m_.clear();
  for (uint k = 0; k < rhs_.size(); k++)
    rhs_[k] = 0.0;

  // update coupling field
  getField("u")->setValues(u_);
}

void ProblemFD1D::assemblyHeat()
{
  for (uint k = 0U; k < n_; k++)
  {
    uint const kLeft = (k != 0) ? k - 1 : k + 1;
    uint const kRight = (k != n_ - 1) ? k + 1 : k - 1;
    m_.add(k, k, 1. / dt_ + alpha_ * 2. / (h_ * h_));
    m_.add(k, kLeft, -alpha_ / (h_ * h_));
    m_.add(k, kRight, -alpha_ / (h_ * h_));
    rhs_[k] = uOld_[k] / dt_ + q_[k];
  }
}

void ProblemFD1D::assemblyHeatCoupled()
{
  // std::vector<double> uExt(n_, 2.0);
  std::vector<double> uExt(n_);
  double const * dataPtr = getField(nameExt_)->dataPtr();
  std::copy(dataPtr, dataPtr + n_, uExt.data());

  double const kAmpli = 10.;
  for (uint k = 1U; k < n_ - 1; k++)
  {
    // matrix
    m_.diag()[k] = 1. / dt_                  // time
                   + alpha_ * 2. / (h_ * h_) // diffusion
                   + kAmpli                  // feedback control
        ;
    m_.diagUp()[k] = -alpha_ / (h_ * h_);   // diffusion
    m_.diagDown()[k] = -alpha_ / (h_ * h_); // diffusion

    // rhs
    rhs_[k] = uOld_[k] / dt_     // time
              + q_[k]            // source
              + kAmpli * uExt[k] // feedback control
        ;
  }
}

// TODO: unify with MatrixCRS version using templates
double computeResidual(
    MatrixTriDiag const & m,
    std::vector<double> const & x,
    std::vector<double> const & b,
    double const area)
{
  std::vector<double> res(x.size());
  auto const tmp = m * x;
  for (uint k = 0; k < res.size(); k++)
    res[k] = (b[k] - tmp[k]);
  double const resNorm = std::sqrt(norm2sq(res) * area);
  return resNorm;
}

std::pair<uint, double> ProblemFD1D::solveTriDiag()
{
  std::vector<double> upPrime(n_);
  std::vector<double> rhsPrime(n_);

  upPrime[0] = m_.diagUp()[0] / m_.diag()[0];
  for (uint k = 1U; k < n_ - 1; k++)
  {
    upPrime[k] = m_.diagUp()[k] / (m_.diag()[k] - m_.diagDown()[k] * upPrime[k - 1]);
  }

  rhsPrime[0] = rhs_[0] / m_.diag()[0];
  for (uint k = 1U; k < n_; k++)
  {
    rhsPrime[k] = (rhs_[k] - m_.diagDown()[k] * rhsPrime[k - 1]) /
                  (m_.diag()[k] - m_.diagDown()[k] * upPrime[k - 1]);
  }

  u_[n_ - 1] = rhsPrime[n_ - 1];
  for (int k = n_ - 2; k >= 0; k--)
  {
    u_[k] = rhsPrime[k] - upPrime[k] * u_[k + 1];
  }

  return std::pair{1U, computeResidual(m_, u_, rhs_, h_)};
}

std::pair<uint, double> ProblemFD1D::solveVanka()
{
  uint maxIter = 1000U;
  double const toll = 1e-6;

  for (uint n = 0U; n < maxIter; n++)
  {
    double const resNorm = computeResidual(m_, u_, rhs_, h_);

    // fmt::print("iter: {:3d}, current residual: {:.8e}\n", j, resNorm);
    if (resNorm < toll)
    {
      return std::pair{n, resNorm};
    }

    // for (uint k = 0U; k < n_; k++)
    //   uOld_[k] = u_[k];

    for (uint k = 1U; k < n_ - 1; k += 2U)
    {
      u_[k] = (rhs_[k] - m_.diagDown()[k] * u_[k - 1] - m_.diagUp()[k] * u_[k + 1]) /
              m_.diag()[k];
    }
    for (uint k = 2U; k < n_ - 1; k += 2U)
    {
      u_[k] = (rhs_[k] - m_.diagDown()[k] * u_[k - 1] - m_.diagUp()[k] * u_[k + 1]) /
              m_.diag()[k];
    }
    u_[0] = (rhs_[0] - m_.diagUp()[0] * u_[1]) / m_.diag()[0];
    u_[n_ - 1] =
        (rhs_[n_ - 1] - m_.diagDown()[n_ - 1] * u_[n_ - 2]) / m_.diag()[n_ - 1];
  }

  return std::pair{maxIter, computeResidual(m_, u_, rhs_, h_)};
}

void ProblemFD1D::print()
{
  std::FILE * out =
      std::fopen((outFile_ + "." + std::to_string(it) + ".dat").c_str(), "w");
  for (uint k = 0; k < u_.size(); k++)
  {
    fmt::print(out, "{:.6e} {:.6e}\n", start_ + k * h_, u_[k]);
  }
  std::fclose(out);

  getField("u")->printVTK(time, it);
}

std::unordered_map<EQN_TYPE, ProblemFD1D::Assembly_T> ProblemFD1D::assemblies_ = {
    {EQN_TYPE::HEAT, [](ProblemFD1D * p) { p->assemblyHeat(); }},
    {EQN_TYPE::HEAT_COUPLED, [](ProblemFD1D * p) { p->assemblyHeatCoupled(); }},
};

std::unordered_map<FD_SOLVER_TYPE, ProblemFD1D::Solver_T> ProblemFD1D::solvers_ = {
    {FD_SOLVER_TYPE::TRIDIAG, [](ProblemFD1D * p) { return p->solveTriDiag(); }},
    {FD_SOLVER_TYPE::VANKA1D, [](ProblemFD1D * p) { return p->solveVanka(); }},
};
