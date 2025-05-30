#pragma once

// std
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <variant>

// local
#include "coupling/field_coupling.hpp"
#include "coupling/mesh_coupling.hpp"
#include "enums.hpp"

namespace cocoa
{

struct Problem
{
  using ConfigList_T = std::unordered_map<
      std::string,
      std::variant<int, std::string, std::vector<std::string>, std::filesystem::path>>;

  Problem() = default;
  Problem(PROBLEM_TYPE type, COUPLING_TYPE couplingType)
      : type_{type}
      , couplingType_{couplingType}
  {}
  virtual ~Problem() = default;

  virtual void setup(ConfigList_T const & configs) = 0;
  virtual bool run() const = 0;
  virtual void advance() = 0;
  virtual uint solve() = 0;
  virtual void print() = 0;

  virtual FieldCoupling * getField(std::string_view name)
  {
    return fieldsCoupling_.at(std::string{name}).get();
  }

  virtual void setField(std::string_view name, FieldCoupling * field)
  {
    fieldsCoupling_.emplace(std::pair{std::string{name}, field});
  }

  MeshCoupling * getMesh() { return meshCoupling_.get(); }

  static std::unique_ptr<Problem> build(PROBLEM_TYPE problemType, EQN_TYPE eqnType);
  static std::unique_ptr<Problem>
  build(std::string_view problemType, std::string_view eqnType);

  PROBLEM_TYPE type_ = PROBLEM_TYPE::NONE;
  COUPLING_TYPE couplingType_ = COUPLING_TYPE::NONE;
  double time = 0.0;
  uint it = 0u;
  uint printStep_ = 1u;
  std::unique_ptr<MeshCoupling> meshCoupling_;
  std::unordered_map<std::string, std::unique_ptr<MeshCoupling>> bcMeshCoupling_;
  std::unordered_map<std::string, std::unique_ptr<FieldCoupling>> fieldsCoupling_;
  bool debug_ = false;
};

} // namespace cocoa
