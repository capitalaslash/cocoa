#pragma once

#include "plugins.hpp"

#ifdef COCOA_ENABLE_OFORG

// fmt
#include <fmt/core.h>

// local
#include "field_coupling.hpp"
#include "problem.hpp"

struct ProblemOForg: public Problem
{
  using ParamList_T = Problem::ParamList_T;

  ProblemOForg() = default;
  ~ProblemOForg() = default;

  void setup(ParamList_T const & params) override
  {
    fmt::print("ProblemOForg::setup() -> init mesh\n");
  }
  bool run() override { return false; }
  void advance() override {}
  void solve() override {}
  void print() override {}
  FieldCoupling getField(std::string_view name) override { return FieldCoupling{}; }
  void setField(std::string_view name, FieldCoupling const & field) override {}
};

#endif

