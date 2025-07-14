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

namespace cocoa
{

struct ProblemFD2D: public Problem
{
  using Assembly_T = std::function<void(ProblemFD2D *)>;
  using Matrix_T = MatrixCSR;

  ProblemFD2D();
  virtual ~ProblemFD2D();

  void setup(Problem::ConfigList_T const & configs) override;
  bool run() const override;
  void advance() override;
  uint solve() override;
  void print() override;
  void printFields();

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
  void assemblyHeatOC();
  void assemblyHeatCoupled();

  void printSetup(std::string_view filename);

  std::string name_;
  MeshFD2D mesh_;
  uint nVars_;
  std::vector<std::string> varNames_;
  VectorFD u_;
  VectorFD uOld_;
  std::unordered_map<std::string, VectorFD> fields_;
  VectorFD q_;
  bool computeCFL_ = false;
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
  std::filesystem::path outputPrefix_ = "./output_fd2d/";
  std::string nameExt_ = "uExternal";
  std::unordered_map<EQN_TYPE, Assembly_T> assemblies_;
  std::unique_ptr<Assembly_T> preSolveFun_ = nullptr;

  static std::unordered_map<FD_SOLVER_TYPE, Solver_T<Matrix_T>> solvers_;
};

const std::vector<uint> sideDOF(std::array<uint, 2U> const & n, FD_BC_SIDE const s);

static constexpr inline uint
sideOffset(std::array<uint, 2U> const & n, FD_BC_SIDE const s)
{
  switch (s)
  {
    using enum FD_BC_SIDE;
  case BOTTOM:
    return +n[0];
  case RIGHT:
    return -1;
  case TOP:
    return -n[0];
  case LEFT:
    return +1;
  default:
    std::abort();
  }
  return 0u;
}

} // namespace cocoa
