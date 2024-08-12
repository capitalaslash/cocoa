#pragma once

#include <fmt/core.h>

#include "field_coupling.hpp"
#include "problem.hpp"

struct ProblemFEMUS: Problem
{
  using ParamList_T = Problem::ParamList_T;

  void setup(ParamList_T const & params) override
  {
    fmt::print("ProblemFEMUS::setup() -> init mesh\n");
  }
  bool run() override { return false; }
  void solve() override {}
  virtual FieldCoupling getField(std::string_view name) override { return FieldCoupling{}; }
  void setField(std::string_view name, FieldCoupling const & field) override {}

  // FEMUSMesh mesh_;
  // std::vector<FEMUSField> fields_;
};

