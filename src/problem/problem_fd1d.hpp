#pragma once

// std
#include <functional>
#include <string>
#include <type_traits>
#include <unordered_map>

// local
#include "enums.hpp"
#include "la.hpp"
#include "problem/fdutils.hpp"
#include "problem/problem.hpp"

struct ProblemFD1D: public Problem
{
  using Assembly_T = std::function<void(ProblemFD1D *)>;
  using Matrix_T = MatrixTriDiag;

  ProblemFD1D();
  virtual ~ProblemFD1D();

  void setup(Problem::ConfigList_T const & configs) override;
  bool run() override;
  void advance() override;
  uint solve() override;
  void print() override;

  void initMeshCoupling();
  void initFieldCoupling();

  void assemblyHeat();
  void assemblyHeatCoupled();

  std::string name_;
  double start_;
  double h_;
  uint n_;
  std::string varName_;
  VectorFD u_;
  VectorFD uOld_;
  VectorFD q_;
  double alpha_;
  double finalTime_;
  double dt_;
  Matrix_T m_;
  VectorFD rhs_;
  FD_SOLVER_TYPE solverType_ = FD_SOLVER_TYPE::TRIDIAG;
  EQN_TYPE eqnType_ = EQN_TYPE::NONE;
  FDBC bcStart_;
  FDBC bcEnd_;
  bool cleanOutput_ = false;
  std::filesystem::path outputPrefix_ = "./output_fd1d";
  std::string nameExt_ = "uExternal";
  std::unordered_map<EQN_TYPE, Assembly_T> assemblies_;

  static std::unordered_map<FD_SOLVER_TYPE, Solver_T<Matrix_T>> solvers_;
};
