#pragma once

#include <fmt/core.h>

#include "med_field.hpp"
#include "problem.hpp"

struct FEMUSMesh
{};

struct FEMUSField
{};

struct ProblemFEMUS: Problem
{
  using ParamList_T = Problem::ParamList_T;

  void setup(ParamList_T const & params) override
  {
    fmt::print("ProblemFEMUS::setup() -> init mesh\n");
  }
  bool run() override { return false; }
  void solve() override {}
  virtual MEDField get_field(std::string_view name) override { return MEDField{}; }
  void set_field(std::string_view name, MEDField const & field) override {}

  FEMUSMesh mesh_;
  std::vector<FEMUSField> fields_;
};
