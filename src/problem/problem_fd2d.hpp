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

struct ProblemFD2D: public Problem
{
  using Vec2D_T = std::array<double, 2U>;
  using Assembly_T = std::function<void(ProblemFD2D *)>;
  using Solver_T = std::function<void(ProblemFD2D *)>;

  ProblemFD2D(): Problem{PROBLEM_TYPE::FD1D, COUPLING_TYPE::NONE} {}
  virtual ~ProblemFD2D() = default;

  void setup(Problem::ParamList_T const & params) override;
  bool run() override;
  void advance() override;
  void solve() override;
  void print() override;

  void initMeshCoupling(std::filesystem::path const & fileName);
  void initFieldCoupling(std::filesystem::path const & fileName);
  void assemblyHeat();
  void assemblyHeatCoupled();
  void solveVankaSCI();
  void solveVankaCB();

  std::string name_;
  Vec2D_T start_;
  Vec2D_T h_;
  std::array<uint, 2U> n_;
  std::vector<double> u_;
  std::vector<double> uOld_;
  std::vector<double> q_;
  double alpha_;
  double finalTime_;
  double dt_;
  MatrixCSR m_;
  std::vector<double> rhs_;
  FD_SOLVER_TYPE solverType_ = FD_SOLVER_TYPE::VANKA2DCB;
  EQN_TYPE eqnType_ = EQN_TYPE::NONE;
  std::array<FDBC, 4U> bcs_;
  std::string outFile_ = "./fd2d";
  std::string nameExt_ = "uExternal";

  static std::unordered_map<EQN_TYPE, Assembly_T> assemblies_;
  static std::unordered_map<FD_SOLVER_TYPE, Solver_T> solvers_;
};
