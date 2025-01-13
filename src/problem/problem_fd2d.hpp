#pragma once

// std
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

// local
#include "enums.hpp"
#include "la.hpp"
#include "problem/fdutils.hpp"
#include "problem/problem.hpp"

struct ProblemFD2D: public Problem
{
  using Vec2D_T = std::array<double, 2U>;
  using Assembly_T = std::function<void(ProblemFD2D *)>;
  using Matrix_T = MatrixCSR;

  ProblemFD2D();
  virtual ~ProblemFD2D();

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
  Vec2D_T start_;
  Vec2D_T h_;
  std::array<uint, 2U> n_;
  std::string varName_;
  std::vector<double> u_;
  std::vector<double> uOld_;
  std::vector<double> q_;
  double alpha_;
  std::array<std::vector<double>, 2> c_;
  double finalTime_;
  double dt_;
  Matrix_T m_;
  std::vector<double> rhs_;
  FD_SOLVER_TYPE solverType_ = FD_SOLVER_TYPE::GAUSSSEIDEL;
  double toll_ = 1.e-6;
  uint maxIters_;
  EQN_TYPE eqnType_ = EQN_TYPE::NONE;
  std::array<FDBC, 4U> bcs_;
  std::filesystem::path outputPrefix_ = "./output_fd2d/";
  std::string nameExt_ = "uExternal";
  std::unordered_map<EQN_TYPE, Assembly_T> assemblies_;

  static std::unordered_map<FD_SOLVER_TYPE, Solver_T<Matrix_T>> solvers_;
};

const std::vector<uint> sideDOF(std::array<uint, 2U> const & n, uint const side);
