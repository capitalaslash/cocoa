#pragma once

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

//- Set the initial time-step according to the solver maxDeltaT
void setDeltaT(Foam::Time & runTime, const Foam::solver & solver);

//- Adjust the time-step according to the solver maxDeltaT
void adjustDeltaT(Foam::Time & runTime, const Foam::solver & solver);

enum struct OFFIELD_TYPE : uint8_t
{
  NONE = 0,
  SCALAR,
  VECTOR,
};

struct ProblemOForg: public Problem
{
  ProblemOForg(): Problem{PROBLEM_TYPE::OFORG, COUPLING_TYPE::NONE} {}
  ~ProblemOForg() { Foam::Info << "End\n" << Foam::endl; }

  void setup(Problem::ConfigList_T const & configs) override;
  bool run() override;
  void advance() override;
  uint solve() override;
  void print() override;

  void initMeshMED(std::string_view name);
  void initFieldMED(std::string_view name, std::filesystem::path path);
  template <typename Field>
  void setDataMED(std::string_view name, Field const & field);

  std::filesystem::path prefix_;
  std::unique_ptr<Foam::Time> runTime_;
  Foam::word solverName_;
  std::unique_ptr<Foam::fvMesh> mesh_;
  std::unique_ptr<Foam::volScalarField> field_;
  Foam::autoPtr<Foam::solver> solverPtr_;
  std::unique_ptr<Foam::pimpleSingleRegionControl> pimple_;
  std::vector<std::pair<std::string, OFFIELD_TYPE>> namesExport_;
  double lastPrint_ = 0.0;
};

#endif
