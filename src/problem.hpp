#pragma once

// std
#include <filesystem>
#include <string>
#include <unordered_map>

// local
#include "field_coupling.hpp"
#include "mesh_coupling.hpp"

struct Problem
{
  using ParamList_T = std::unordered_map<std::string, std::filesystem::path>;

  Problem() = default;
  explicit Problem(PROBLEM_TYPE type): type_{type} {}
  virtual ~Problem() = default;

  virtual void setup(ParamList_T const & params) = 0;
  virtual bool run() = 0;
  virtual void advance() = 0;
  virtual void solve() = 0;
  virtual void print() = 0;

  virtual FieldCoupling * getField(std::string_view name)
  {
    return fieldsCoupling_.at(std::string{name}).get();
  }

  virtual void setField(std::string_view name, FieldCoupling * field)
  {
    fieldsCoupling_.at(std::string{name}).reset(field);
  }

  MeshCoupling * getMesh() { return meshCoupling_.get(); }

  PROBLEM_TYPE type_ = PROBLEM_TYPE::NONE;
  COUPLING_TYPE couplingType_ = COUPLING_TYPE::NONE;
  double time = 0.0;
  int it = 0;
  std::unique_ptr<MeshCoupling> meshCoupling_;
  std::unordered_map<std::string, std::unique_ptr<FieldCoupling>> fieldsCoupling_;
};
