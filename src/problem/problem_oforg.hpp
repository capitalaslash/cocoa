#pragma once

#include "coupling/field_coupling.hpp"
#include "coupling/mesh_coupling.hpp"
#include "plugins.hpp"

#ifdef COCOA_ENABLE_OFORG

// fmt
#include <fmt/core.h>

// openfoam.org
#include <argList.H>
#include <pimpleSingleRegionControl.H>
#include <solver.H>

// local
#include "problem/problem.hpp"

namespace cocoa
{

//- Set the initial time-step according to the solver maxDeltaT
void setDeltaT(Foam::Time & runTime, const Foam::solver & solver);

//- Adjust the time-step according to the solver maxDeltaT
void adjustDeltaT(Foam::Time & runTime, const Foam::solver & solver);

enum struct OFFIELD_TYPE : uint8_t
{
  NONE = 0u,
  SCALAR,
  VECTOR,
};

struct ProblemOForg: public Problem
{
  ProblemOForg(): Problem{PROBLEM_TYPE::OFORG} {}
  ~ProblemOForg() { Foam::Info << "End\n" << Foam::endl; }

  void setup(Problem::ConfigList_T const & configs) override;
  bool run() const override;
  void advance() override;
  uint solve() override;
  void print() override;
  void printVTK();

  // TODO: get variables from mesh
  std::vector<std::string> varNames() override { return {"p", "U"}; }
  Marker findRegion(std::string_view name) override;

  std::unique_ptr<MeshCoupling> initMeshCoupling(COUPLING_TYPE type) override;
  std::unique_ptr<FieldCoupling> initFieldCoupling(
      COUPLING_TYPE type, std::string_view name, MeshCoupling const * mesh) override;
  void setFieldData(FieldCoupling * field) override;
  void getFieldData(FieldCoupling const & field) override;

  std::filesystem::path prefix_;
  std::unique_ptr<Foam::Time> runTime_;
  Foam::word solverName_;
  std::unique_ptr<Foam::fvMesh> mesh_;
  // std::unique_ptr<Foam::volScalarField> field_;
  Foam::autoPtr<Foam::solver> solverPtr_;
  std::unique_ptr<Foam::pimpleSingleRegionControl> pimple_;
  std::filesystem::path outputVTK_;
  std::unique_ptr<MeshCoupling> meshExport_;
  std::unordered_map<std::string, std::unique_ptr<FieldCoupling>> fieldsExport_;
  double lastPrint_ = 0.0;

  static bool argInit;
};

} // namespace cocoa

#endif
