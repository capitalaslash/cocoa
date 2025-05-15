#pragma once

#include "plugins.hpp"

#ifdef COCOA_ENABLE_FEMUS

// fmt
#include <fmt/core.h>

// local
#include "field_coupling.hpp"
#include "problem.hpp"

namespace cocoa
{

struct ProblemFEMUS: public Problem
{
  ProblemFEMUS: Problem{PROBLEM_TYPE::FEMUS} {}
  ~ProblemFEMUS() = default;

  void setup(Probelm::ConfigList_T const & configs) override
  {
    fmt::println("ProblemFEMUS::setup() -> init mesh");
  }
  bool run() override { return false; }
  uint solve() override {}

  // FEMUSMesh mesh_;
  // std::vector<FEMUSField> fields_;
};

} // namespace cocoa

#endif
