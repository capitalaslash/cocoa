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
  void assemblyHeatOC();
  void assemblyHeatCoupled();

  void printSetup(std::string_view filename);

  std::string name_;
  MeshFD2D mesh_;
  uint nVars_;
  std::vector<std::string> varNames_;
  VectorFD u_;
  VectorFD uOld_;
  VectorFD q_;
  bool computeCFL_ = false;
  std::array<VectorFD, 2u> c_;
  ParamsFD params_;
  double finalTime_;
  double dt_;
  Matrix_T m_;
  VectorFD rhs_;
  FD_SOLVER_TYPE solverType_ = FD_SOLVER_TYPE::GAUSS_SEIDEL;
  double tol_ = 1.e-6;
  uint maxIters_;
  EQN_TYPE eqnType_ = EQN_TYPE::NONE;
  std::vector<FDBCList2D> bcs_;
  bool cleanOutput_ = false;
  uint printStep_ = 1u;
  std::filesystem::path outputPrefix_ = "./output_fd2d/";
  std::string nameExt_ = "uExternal";
  std::unordered_map<EQN_TYPE, Assembly_T> assemblies_;

  static std::unordered_map<FD_SOLVER_TYPE, Solver_T<Matrix_T>> solvers_;
};

const std::vector<uint> sideDOF(std::array<uint, 2U> const & n, FD_BC_SIDE const s);
