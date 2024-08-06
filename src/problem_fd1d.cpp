#include "problem_fd1d.hpp"

#include <filesystem>
#include <fstream>

#include <fmt/core.h>

void ProblemFD1D::setup(Problem::ParamList_T const & params)
{
  // default values
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
  // bcs
  bcStart_ = 1.0;
  bcEnd_ = 0.0;

  // read configuration from file
  std::filesystem::path configFile = params.at("config_file");
  std::ifstream in(configFile, std::ios::in);
  std::string buffer;
  while (std::getline(in, buffer, '\n'))
  {
    std::istringstream bufferStream{buffer};
    std::string token;
    while (std::getline(bufferStream, token, ' '))
    {
      if (token == "start:")
        bufferStream >> start;
      else if (token == "end:")
        bufferStream >> end;
      else if (token == "n_elems:")
        bufferStream >> nElems;
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
      else if (token == "bc_start:")
        bufferStream >> bcStart_;
      else if (token == "bc_end:")
        bufferStream >> bcEnd_;
      else if (token == "out_file:")
        bufferStream >> outFile_;
      else
      {
        fmt::print(stderr, "key {} invalid\n", token);
        bufferStream >> token;
      }
    }
  }

  // mesh
  auto const h = (end - start) / nElems;
  points_.resize(nElems + 1);
  for (uint k = 0; k < points_.size(); k++)
  {
    points_[k] = start + k * h;
  }
  auto const dirPath = std::filesystem::path(outFile_).parent_path();
  initMeshMED(dirPath / "mesh");

  // fields
  u_.resize(points_.size(), uInit);
  uOld_.resize(points_.size(), uInit);
  q_.resize(points_.size(), qValue);
  initFieldMED(dirPath / "u");

  // io
  std::filesystem::create_directories(dirPath);
}

void ProblemFD1D::initMeshMED(std::filesystem::path const & fileName)
{
  // coords format: x_0, y_0, z_0, x_1, ...
  std::vector<double> coords(points_.size() * 3);
  for (uint k = 0; k < points_.size(); k++)
  {
    coords[3 * k] = points_[k];
    coords[3 * k + 1] = 0.0;
    coords[3 * k + 2] = 0.0;
  }

  // conn format: elem0_numpts, id_0, id_1, ..., elem1_numpts, ...
  auto const nElems = points_.size() - 1;
  std::vector<mcIdType> conn(nElems * 3);
  for (uint k = 0; k < nElems; k++)
  {
    conn[3 * k] = MEDCellTypeToIKCell(MED_CELL_TYPE::LINE2);
    conn[3 * k + 1] = k;
    conn[3 * k + 2] = k + 1;
  }

  // offsets format: sum_0^k elemk_numpts + 1,
  std::vector<mcIdType> offsets(nElems + 1);
  offsets[0] = 0;
  for (uint k = 0; k < nElems; k++)
  {
    offsets[k + 1] = offsets[k] + 3;
  }

  meshMED_.init(fileName.string(), coords, conn, offsets);
}

void ProblemFD1D::initFieldMED(std::filesystem::path const & fileName)
{
  uMED_.init(fileName.filename().string(), meshMED_);
  uMED_.setValues(u_);
  uMED_.initIO(fileName.string() + "_med.");
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
  fmt::print("fd1d::solve(), time = {:.6e}, dt = {:.6e}\n", time, dt_);
  uint const n = points_.size();
  double const h = 1. / (n - 1);

  // update
  for (uint k = 0; k < n; k++)
    uOld_[k] = u_[k];

  Matrix m{n};
  std::vector<double> rhs(n);

  // bc start
  m.diag[0] = 1.0;
  m.diagUp[0] = 0.0;
  rhs[0] = bcStart_;

  // bc end
  m.diag[n - 1] = 1.0;
  m.diagDown[n - 1] = 0.0;
  rhs[n - 1] = bcEnd_;

  // assembly
  for (uint k = 1U; k < n - 1; k++)
  {
    m.diag[k] = 1. / dt_ + diff_ * 2 / (h * h);
    m.diagUp[k] = -diff_ / (h * h);
    m.diagDown[k] = -diff_ / (h * h);
    rhs[k] = q_[k] + uOld_[k] / dt_;
  }

  // solve
  std::vector<double> upPrime(n);
  std::vector<double> rhsPrime(n);

  upPrime[0] = m.diagUp[0] / m.diag[0];
  for (uint k = 1U; k < n - 1; k++)
  {
    upPrime[k] = m.diagUp[k] / (m.diag[k] - m.diagDown[k] * upPrime[k - 1]);
  }

  rhsPrime[0] = rhs[0] / m.diag[0];
  for (uint k = 1U; k < n; k++)
  {
    rhsPrime[k] = (rhs[k] - m.diagDown[k] * rhsPrime[k - 1]) /
                  (m.diag[k] - m.diagDown[k] * upPrime[k - 1]);
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
    fmt::print(out, "{:.6e} {:.6e}\n", points_[k], u_[k]);
  }
  std::fclose(out);

  uMED_.setValues(u_);
  uMED_.printVTK(time, it);
}
