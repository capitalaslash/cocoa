#pragma once

#include "plugins.hpp"

#ifdef COCOA_ENABLE_FEMUS

// fmt
#include <fmt/core.h>

// local
#include "field_coupling.hpp"
#include "problem.hpp"

struct ProblemFEMUS: public Problem
{
  ProblemFEMUS: Problem{PROBLEM_TYPE::FEMUS} {}
  ~ProblemFEMUS() = default;

  using ParamList_T = Problem::ParamList_T;

  void setup(ParamList_T const & params) override
  {
    fmt::print("ProblemFEMUS::setup() -> init mesh\n");
  }
  bool run() override { return false; }
  void solve() override {}

  // FEMUSMesh mesh_;
  // std::vector<FEMUSField> fields_;
};

#endif
