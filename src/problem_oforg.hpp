#pragma once

#include <fmt/core.h>

#include "field_coupling.hpp"
#include "problem.hpp"

struct ProblemOForg: Problem
{
  using ParamList_T = Problem::ParamList_T;

  void setup(ParamList_T const & params) override
  {
    fmt::print("ProblemOForg::setup() -> init mesh\n");
  }
  bool run() override { return false; }
  void solve() override {}
  FieldCoupling getField(std::string_view name) override { return FieldCoupling{}; }
  void setField(std::string_view name, FieldCoupling const & field) override {}

  // OForgMesh mesh_;
  // std::vector<OForgField> fields_;
};
