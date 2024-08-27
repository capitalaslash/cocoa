#include "problem_fd1d.hpp"

// std
#include <cassert>
#include <filesystem>
#include <fstream>
#include <unordered_map>

// libfmt
#include <fmt/core.h>

// local
#include "field_factory.hpp"
#include "mesh_med.hpp"

void ProblemFD1D::setup(Problem::ParamList_T const & params)
{
  // default values
  name_ = "empty";
  // mesh
  double start = 0.0;
  double end = 0.0;
  uint nElems = 10U;
  // fields
  double uInit = 0.0;
  double qValue = 1.0;
  // physical constants
  diff_ = 0.2;
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
      else if (token == "diff:")
        bufferStream >> diff_;
      else if (token == "start_time:")
        bufferStream >> time;
      else if (token == "final_time:")
        bufferStream >> finalTime_;
      else if (token == "dt:")
        bufferStream >> dt_;
      else if (token == "assembly_name:")
        bufferStream >> assemblyName_;
      else if (token == "bc_start_type:")
      {
        bufferStream >> token;
        bcStartType_ = str2fdbc(token);
      }
      else if (token == "bc_start_value:")
        bufferStream >> bcStartValue_;
      else if (token == "bc_end_type:")
      {
        bufferStream >> token;
        bcEndType_ = str2fdbc(token);
      }
      else if (token == "bc_end_value:")
        bufferStream >> bcEndValue_;
      else if (token == "out_file:")
        bufferStream >> outFile_;
      else
      {
        fmt::print(stderr, "key {} invalid\n", token);
        bufferStream >> token;
      }
    }
    if (assemblyName_ != "")
    {
      assert(assemblies_.contains(assemblyName_));
    }
  }

  // mesh
  start_ = start;
  h_ = (end - start) / nElems;
  n_ = nElems + 1;
  auto const dirPath = std::filesystem::path(outFile_).parent_path();
  initMeshMED(dirPath / "mesh");

  // fields
  u_.resize(n_, uInit);
  uOld_.resize(n_, uInit);
  q_.resize(n_, qValue);
  initFieldMED(dirPath / "u");

  // linear algebra
  m_.init(n_);
  rhs_.resize(n_);

  // io
  std::filesystem::create_directories(dirPath);
}

void ProblemFD1D::initMeshMED(std::filesystem::path const & fileName)
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
    conn[3 * k] = 2; // MEDCellTypeToIKCell(MED_CELL_TYPE::LINE2);
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

  // meshCoupling_.reset(buildMesh(couplingType_));
  meshCoupling_.reset(new MeshMED);
  meshCoupling_->init(fileName.string(), 1U, 1U, coords, conn, offsets);
}

void ProblemFD1D::initFieldMED(std::filesystem::path const & fileName)
{
  auto [kvPairU, successU] = fieldsCoupling_.emplace("u", buildField(couplingType_));
  assert(successU);
  kvPairU->second->init(fileName.filename().string(), meshCoupling_.get());
  kvPairU->second->setValues(u_);
  kvPairU->second->initIO(fileName.string() + "_med.");

  auto [kvPairExt, successExt] =
      fieldsCoupling_.emplace(nameExt_, buildField(couplingType_));
  assert(successExt);
  kvPairExt->second->init(nameExt_, meshCoupling_.get());
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
  fmt::print("{}, time = {:.6e}, dt = {:.6e}\n", name_, time, dt_);
  uint const n = n_;
  double const h = 1. / (n - 1);

  // update
  for (uint k = 0; k < n; k++)
    uOld_[k] = u_[k];

  // bc start
  switch (bcStartType_)
  {
  case FDBC_TYPE::DIRICHLET:
  {
    m_.diag[0] = 1.0;
    m_.diagUp[0] = 0.0;
    rhs_[0] = bcStartValue_;
    break;
  }
  case FDBC_TYPE::NEUMANN:
  {
    m_.diag[0] = 1.0;
    m_.diagUp[0] = -1.0;
    rhs_[0] = bcStartValue_ * h;
    break;
  }
  default:
  {
    fmt::print(stderr, "no bc start specified!\n");
    std::abort();
  }
  }

  // bc end
  switch (bcEndType_)
  {
  case FDBC_TYPE::DIRICHLET:
  {
    m_.diag[n - 1] = 1.0;
    m_.diagDown[n - 1] = 0.0;
    rhs_[n - 1] = bcEndValue_;
    break;
  }
  case FDBC_TYPE::NEUMANN:
  {
    m_.diag[n - 1] = -1.0;
    m_.diagDown[n - 1] = 1.0;
    rhs_[n - 1] = bcEndValue_ * h;
    break;
  }
  default:
  {
    fmt::print(stderr, "no bc end specified!\n");
    std::abort();
  }
  }

  // assembly
  // std::vector<double> uExt(n);
  // double const * dataPtr =
  //     getField(nameExt_)->fieldPtr_->getArray()->getConstPointer();
  // std::copy(dataPtr, dataPtr + n, uExt.data());
  // double const kAmpli = 0.1;
  // double const h = 1. / (n - 1);
  // for (uint k = 1U; k < n - 1; k++)
  // {
  //   m.diag[k] = 1. / dt_ + diff_ * 2. / (h * h);
  //   m.diagUp[k] = -diff_ / (h * h);
  //   m.diagDown[k] = -diff_ / (h * h);
  //   rhs[k] = uOld_[k] / dt_ + q_[k] + couple_ * kAmpli * (uOld_[k] - uExt[k]);
  // }
  assemblies_.at(assemblyName_)(this);

  // solve
  std::vector<double> upPrime(n);
  std::vector<double> rhsPrime(n);

  upPrime[0] = m_.diagUp[0] / m_.diag[0];
  for (uint k = 1U; k < n - 1; k++)
  {
    upPrime[k] = m_.diagUp[k] / (m_.diag[k] - m_.diagDown[k] * upPrime[k - 1]);
  }

  rhsPrime[0] = rhs_[0] / m_.diag[0];
  for (uint k = 1U; k < n; k++)
  {
    rhsPrime[k] = (rhs_[k] - m_.diagDown[k] * rhsPrime[k - 1]) /
                  (m_.diag[k] - m_.diagDown[k] * upPrime[k - 1]);
  }

  u_[n - 1] = rhsPrime[n - 1];
  for (int k = n - 2; k >= 0; k--)
  {
    u_[k] = rhsPrime[k] - upPrime[k] * u_[k + 1];
  }
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

  auto uCoupling = getField("u");
  uCoupling->setValues(u_);
  uCoupling->printVTK(time, it);
}

std::unordered_map<std::string, ProblemFD1D::Assembly_T> ProblemFD1D::assemblies_;

void setFD1DAssemblies()
{
  ProblemFD1D::assemblies_.emplace(
      "heat",
      [](ProblemFD1D * p)
      {
        for (uint k = 1U; k < p->n_ - 1; k++)
        {
          p->m_.diag[k] = 1. / p->dt_ + p->diff_ * 2. / (p->h_ * p->h_);
          p->m_.diagUp[k] = -p->diff_ / (p->h_ * p->h_);
          p->m_.diagDown[k] = -p->diff_ / (p->h_ * p->h_);
          p->rhs_[k] = p->uOld_[k] / p->dt_ + p->q_[k];
        }
      });

  ProblemFD1D::assemblies_.emplace(
      "heatCoupled",
      [](ProblemFD1D * p)
      {
        // std::vector<double> uExt(p->n_, 2.0);
        std::vector<double> uExt(p->n_);
        double const * dataPtr = p->getField(p->nameExt_)->dataPtr();
        std::copy(dataPtr, dataPtr + p->n_, uExt.data());

        double const kAmpli = 10.;
        for (uint k = 1U; k < p->n_ - 1; k++)
        {
          // matrix
          p->m_.diag[k] = 1. / p->dt_                       // time
                          + p->diff_ * 2. / (p->h_ * p->h_) // diffusion
                          + kAmpli                          // feedback control
              ;
          p->m_.diagUp[k] = -p->diff_ / (p->h_ * p->h_);   // diffusion
          p->m_.diagDown[k] = -p->diff_ / (p->h_ * p->h_); // diffusion

          // rhs
          p->rhs_[k] = p->uOld_[k] / p->dt_ // time
                       + p->q_[k]           // source
                       + kAmpli * uExt[k]   // feedback control
              ;
        }
      });
}
