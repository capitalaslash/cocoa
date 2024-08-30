#pragma once

// std
#include <cassert>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// local
#include "field_coupling.hpp"
#include "problem.hpp"

struct ProblemFD1D: public Problem
{
  using Assembly_T = std::function<void(ProblemFD1D *)>;

  struct Matrix
  {
    Matrix() = default;
    explicit Matrix(uint n): diag(n), diagUp(n), diagDown(n) {}
    ~Matrix() = default;

    void init(uint n)
    {
      diag.resize(n);
      diagUp.resize(n);
      diagDown.resize(n);
    }

    std::vector<double> diag;
    std::vector<double> diagUp;
    std::vector<double> diagDown;
  };

  ProblemFD1D(): Problem{PROBLEM_TYPE::FD1D, COUPLING_TYPE::NONE} {}
  ~ProblemFD1D() = default;

  void setup(Problem::ParamList_T const & params) override;
  bool run() override;
  void advance() override;
  void solve() override;
  void print() override;

  void initMeshCoupling(std::filesystem::path const & fileName);
  void initFieldCoupling(std::filesystem::path const & fileName);

  std::string name_;
  double start_;
  double h_;
  uint n_;
  std::vector<double> u_;
  std::vector<double> uOld_;
  std::vector<double> q_;
  double diff_;
  double finalTime_;
  double dt_;
  Matrix m_;
  std::vector<double> rhs_;
  EQN_TYPE eqnType_;
  FDBC_TYPE bcStartType_ = FDBC_TYPE::NONE;
  double bcStartValue_ = 0.0;
  FDBC_TYPE bcEndType_ = FDBC_TYPE::NONE;
  double bcEndValue_ = 0.0;
  std::string outFile_ = "./fd1d";
  std::string nameExt_ = "uExternal";

  static std::unordered_map<EQN_TYPE, Assembly_T> assemblies_;
};
