#pragma once

// std
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

// local
#include "enums.hpp"
#include "problem/fdutils.hpp"
#include "problem/problem.hpp"

struct ProblemFD1D: public Problem
{
  using Assembly_T = std::function<void(ProblemFD1D *)>;
  using Solver_T = std::function<void(ProblemFD1D *)>;

  ProblemFD1D(): Problem{PROBLEM_TYPE::FD1D, COUPLING_TYPE::NONE} {}
  ~ProblemFD1D() = default;

  void setup(Problem::ParamList_T const & params) override;
  bool run() override;
  void advance() override;
  void solve() override;
  void print() override;

  void initMeshCoupling(std::filesystem::path const & fileName);
  void initFieldCoupling(std::filesystem::path const & fileName);
  void assemblyHeat();
  void assemblyHeatCoupled();
  void solveTriDiag();
  void solveVanka();

  std::string name_;
  double start_;
  double h_;
  uint n_;
  std::vector<double> u_;
  std::vector<double> uOld_;
  std::vector<double> q_;
  double alpha_;
  double finalTime_;
  double dt_;
  MatrixTriDiag m_;
  std::vector<double> rhs_;
  FD_SOLVER_TYPE solverType_ = FD_SOLVER_TYPE::TRIDIAG;
  EQN_TYPE eqnType_ = EQN_TYPE::NONE;
  FD_BC_TYPE bcStartType_ = FD_BC_TYPE::NONE;
  double bcStartValue_ = 0.0;
  FD_BC_TYPE bcEndType_ = FD_BC_TYPE::NONE;
  double bcEndValue_ = 0.0;
  std::string outFile_ = "./fd1d";
  std::string nameExt_ = "uExternal";

  static std::unordered_map<EQN_TYPE, Assembly_T> assemblies_;
  static std::unordered_map<FD_SOLVER_TYPE, Solver_T> solvers_;
};