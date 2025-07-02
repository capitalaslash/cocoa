#pragma once

// stl
#include <filesystem>
#include <functional>
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
  Problem(PROBLEM_TYPE type): type_{type} {}
  virtual ~Problem() = default;

  virtual void setup(ConfigList_T const & configs) = 0;
  virtual bool run() const = 0;
  virtual void advance() = 0;
  virtual uint solve() = 0;
  virtual void print() = 0;

  virtual std::vector<std::string> varNames() = 0;
  virtual Marker findRegion(std::string_view name) = 0;

  virtual std::unique_ptr<MeshCoupling> initMeshCoupling(COUPLING_TYPE type) = 0;
  virtual std::unique_ptr<FieldCoupling> initFieldCoupling(
      COUPLING_TYPE type, std::string_view name, MeshCoupling const * mesh) = 0;
  virtual void setFieldData(FieldCoupling * field) = 0;
  virtual void getFieldData(FieldCoupling const & field) = 0;

  static std::unique_ptr<Problem> build(PROBLEM_TYPE problemType, EQN_TYPE eqnType);
  static std::unique_ptr<Problem>
  build(std::string_view problemType, std::string_view eqnType);

  PROBLEM_TYPE type_ = PROBLEM_TYPE::NONE;
  std::unordered_map<std::string, double *> dataPtr_;
  double time = 0.0;
  uint it = 0u;
  uint printStep_ = 1u;
  bool debug_ = false;

  static std::unordered_map<PROBLEM_TYPE, std::function<std::unique_ptr<Problem>(void)>>
      builder_;
};

} // namespace cocoa
