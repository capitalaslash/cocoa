#include "problem/problem_fd1d.hpp"

// std
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <unordered_map>

// libfmt
#include <fmt/core.h>
#include <fmt/ranges.h>

// local
#include "coupling/field_coupling.hpp"
#include "coupling/mesh_coupling.hpp"
#include "enums.hpp"
#include "problem/fdutils.hpp"

void ProblemFD1D::setup(Problem::ConfigList_T const & configs)
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
  std::filesystem::path const configFile = configs.at("config_file");
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
        // this is a comment, consume whole line
        while (bufferStream)
        {
          bufferStream >> token;
        }
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
      else if (token == "field_external:")
        bufferStream >> nameExt_;
      else if (token == "var_name:")
        bufferStream >> varName_;
      else if (token == "initial_value:")
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
      else if (token == "clean_output:")
        bufferStream >> cleanOutput_;
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
  h_ = (end - start) / nElems;
  n_ = nElems + 1;
  initMeshCoupling();

  // fields
  u_.resize(n_, uInit);
  uOld_.resize(n_, uInit);
  q_.resize(n_, qValue);
  initFieldCoupling();

  // linear algebra
  m_.init(n_);
  rhs_.resize(n_);

  // io
  if (cleanOutput_)
    std::filesystem::remove_all(outputPrefix_);
  std::filesystem::create_directories(outputPrefix_);
}

void ProblemFD1D::initMeshCoupling()
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
  meshCoupling_->init("mesh_fd1d", 1U, 1U, coords, conn, offsets);
}

void ProblemFD1D::initFieldCoupling()
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
    m_.clearRow(0);
    m_.add(0, 0, 1.0);
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
    m_.add(n_ - 1, n_ - 1, 1.0);
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

  m_.close();

  // solve
  auto const [numIters, residual] =
      solvers_.at(solverType_)(m_, rhs_, u_, 1.e-6, 1000u);
  fmt::print("num iters: {:4d}, ", numIters);
  double const rhsNorm = std::sqrt(norm2sq(rhs_) * h_);
  fmt::print("relative residual: {:.8e}\n", residual / rhsNorm);

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

void ProblemFD1D::assemblyHeat()
{
  // eqn:
  // du/dt - alpha * d^2u/dx^2 = q
  // discretization:
  // u_m / dt - alpha (u_l - 2 * u_m + u_r) / h^2 =
  // uold_m / dt + q_m
  // grouping:
  // (1 / dt + 2 * alpha / h^2) * u_m
  // - (alpha / h^2) * u_l
  // - (alpha / h^2) * u_r
  // = uold_m / dt + q_m
  for (uint k = 0u; k < n_; k++)
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
  for (uint k = 0u; k < n_; k++)
  {
    uint const kLeft = (k != 0) ? k - 1 : k + 1;
    uint const kRight = (k != n_ - 1) ? k + 1 : k - 1;
    // matrix
    m_.add(
        k,
        k,
        1. / dt_                      // time
            + alpha_ * 2. / (h_ * h_) // diffusion
            + kAmpli);
    m_.add(
        k, kRight, -alpha_ / (h_ * h_) // diffusion
    );
    m_.add(
        k, kLeft, -alpha_ / (h_ * h_) // diffusion
    );

    // rhs
    rhs_[k] = uOld_[k] / dt_     // time
              + q_[k]            // source
              + kAmpli * uExt[k] // feedback control
        ;
  }
}

void ProblemFD1D::print()
{
  auto const filename =
      fmt::format("{}.{}.dat", (outputPrefix_ / varName_).string(), it);
  std::FILE * out = std::fopen(filename.c_str(), "w");
  for (uint k = 0; k < u_.size(); k++)
  {
    fmt::print(out, "{:.6e} {:.6e}\n", start_ + k * h_, u_[k]);
  }
  std::fclose(out);

  getField(varName_)->printVTK(time, it);
}

std::unordered_map<EQN_TYPE, ProblemFD1D::Assembly_T> ProblemFD1D::assemblies_ = {
    {EQN_TYPE::HEAT, [](ProblemFD1D * p) { p->assemblyHeat(); }},
    {EQN_TYPE::HEAT_COUPLED, [](ProblemFD1D * p) { p->assemblyHeatCoupled(); }},
};

std::unordered_map<FD_SOLVER_TYPE, Solver_T<ProblemFD1D::Matrix_T>>
    ProblemFD1D::solvers_ = {
        {FD_SOLVER_TYPE::NONE,
         [](MatrixTriDiag const & m,
            VectorFD const & b,
            VectorFD & x,
            double const residual,
            uint const maxIters) {
           return SolverInfo{0u, 0.0};
         }},
        {FD_SOLVER_TYPE::JACOBI, &solveJacobi<MatrixTriDiag>},
        {FD_SOLVER_TYPE::TRIDIAG,
         [](MatrixTriDiag const & m,
            VectorFD const & b,
            VectorFD & x,
            double const,
            uint const) { return solveTriDiag(m, b, x); }},
        {FD_SOLVER_TYPE::VANKA1D, &solveVanka1D},
};
