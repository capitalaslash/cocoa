#pragma once

#include <fmt/core.h>

#include "med_field.hpp"
#include "problem.hpp"

struct OForgMesh
{};

struct OForgField
{};

struct ProblemOForg: Problem
{
  using ParamList_T = Problem::ParamList_T;

  void setup(ParamList_T const & params) override
  {
    fmt::print("ProblemOForg::setup() -> init mesh\n");
  }
  bool run() override { return false; }
  void solve() override {}
  MEDField getField(std::string_view name) override { return MEDField{}; }
  void setField(std::string_view name, MEDField const & field) override {}

  OForgMesh mesh_;
  std::vector<OForgField> fields_;
};
