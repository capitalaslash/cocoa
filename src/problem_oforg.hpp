#pragma once

#include "plugins.hpp"

#ifdef COCOA_ENABLE_OFORG

// fmt
#include <fmt/core.h>

// local
#include "problem.hpp"

struct ProblemOForg: public Problem
{
  using ParamList_T = Problem::ParamList_T;

  ProblemOForg(): Problem{PROBLEM_TYPE::OFORG, COUPLING_TYPE::NONE} {}
  ~ProblemOForg() = default;

  void setup(ParamList_T const & params) override
  {
    fmt::print("ProblemOForg::setup() -> init mesh\n");
  }
  bool run() override { return false; }
  void advance() override {}
  void solve() override {}
  void print() override {}
};

#endif
