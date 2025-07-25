#pragma once

// std
#include <functional>
#include <string>
#include <string_view>
#include <unordered_map>

// local
#include "enums.hpp"
#include "la.hpp"
#include "problem/fdutils.hpp"
#include "problem/problem.hpp"

namespace cocoa
{

struct ProblemFD1D: public Problem
{
  using Assembly_T = std::function<void(ProblemFD1D *)>;
  using Matrix_T = MatrixCSR;

  ProblemFD1D();
  ~ProblemFD1D();

  void setup(Problem::ConfigList_T const & configs) override;
  bool run() const override;
  void advance() override;
  uint solve() override;
  void print() override;

  std::vector<std::string> varNames() override { return varNames_; }
  Marker findRegion(std::string_view name) override;

  std::unique_ptr<MeshCoupling> initMeshCoupling(
      COUPLING_TYPE type,
      COUPLING_SCOPE scope,
      Marker marker,
      std::string_view bdName) override;
  std::unique_ptr<FieldCoupling> initFieldCoupling(
      COUPLING_TYPE type, std::string_view name, MeshCoupling const * mesh) override;
  void setFieldData(FieldCoupling * field) override;
  void getFieldData(FieldCoupling const & field) override;

  void initOutput();

  void assemblyHeat();
  void assemblyHeatCoupled();
  void assemblyHeatOC();

  std::string name_;
  MeshFD1D mesh_;
  uint nVars_;
  std::vector<std::string> varNames_;
  VectorFD u_;
  VectorFD uOld_;
  std::unordered_map<std::string, VectorFD> fields_;
  ParamsFD params_;
  double finalTime_;
  double dt_;
  Matrix_T m_;
  VectorFD rhs_;
  FD_SOLVER_TYPE solverType_ = FD_SOLVER_TYPE::GAUSS_SEIDEL;
  uint maxIters_ = 1000u;
  double tol_ = 1.0e-6;
  EQN_TYPE eqnType_ = EQN_TYPE::NONE;
  std::vector<FDBCList1D> bcs_;
  bool cleanOutput_ = false;
  std::filesystem::path outputPrefix_ = "./output_fd1d";
  std::unordered_map<EQN_TYPE, Assembly_T> assemblies_;

  static std::unordered_map<FD_SOLVER_TYPE, Solver_T<Matrix_T>> solvers_;
};

} // namespace cocoa
