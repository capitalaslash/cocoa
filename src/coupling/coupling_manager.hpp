#pragma once

// stl
#include <array>
#include <memory>
#include <string_view>
#include <unordered_map>

// local
#include "coupling/field_coupling.hpp"
#include "coupling/mesh_coupling.hpp"
#include "enums.hpp"
#include "problem/problem.hpp"

namespace cocoa
{

struct Problem;

struct CouplingInterface
{
  CouplingInterface() = default;

  CouplingInterface(Problem * p, std::vector<std::string> fn): ptr{p}, fieldNames{fn} {}

  CouplingInterface(Problem * p, std::string_view bdn, std::vector<std::string> fn)
      : ptr{p}
      , marker{p->findRegion(bdn)}
      , bdName(bdn)
      , fieldNames{fn}
  {}

  ~CouplingInterface() = default;

  Problem * ptr = nullptr;
  Marker marker = markerNotSet;
  std::string bdName = "";
  std::vector<std::string> fieldNames;
};

enum class INTERPOLATION_METHOD : uint8_t
{
  NONE = 0u,
  P0P0,
  P1P1,
};

inline constexpr std::string interpolationMethod2str(INTERPOLATION_METHOD method)
{
  using enum INTERPOLATION_METHOD;
  switch (method)
  {
  case P0P0:
    return "P0P0";
  case P1P1:
    return "P1P1";
  default:
    std::abort();
  }
}

struct CouplingManager
{
  CouplingManager() = default;
  CouplingManager(COUPLING_TYPE type, COUPLING_SCOPE scope): type_{type}, scope_{scope}
  {}
  virtual ~CouplingManager() = default;

  virtual void setup(
      CouplingInterface interfaceSrc,
      CouplingInterface interfaceTgt,
      INTERPOLATION_METHOD method) = 0;
  virtual void initFieldCoupling() = 0;

  void getField(std::string_view fieldName)
  {
    auto * field = fields_[0].at(std::string{fieldName}).get();
    interfaces_[0].ptr->setFieldData(field);
  }

  void setField(std::string_view fieldName)
  {
    auto const & field = *(fields_[1].at(std::string{fieldName}));
    interfaces_[1].ptr->getFieldData(field);
  }

  void updateFields()
  {
    for (auto & [name, _]: fields_[0])
    {
      this->getField(name);
    }
  }

  virtual void
  project(std::string_view fieldNameSrc, std::string_view fieldNameTgt) = 0;
  void project(std::string_view fieldName) { this->project(fieldName, fieldName); }

  void printVTK()
  {
    for (auto k = 0u; k < 2u; k++)
    {
      auto p = interfaces_[k].ptr;
      if (p->it % p->printStep_ == 0u)
        for (auto const & [_, field]: fields_[k])
        {
          field->printVTK(p->time, p->it);
        }
    }
  }

  static std::unique_ptr<CouplingManager>
  build(COUPLING_TYPE type, COUPLING_SCOPE scope);
  static std::unique_ptr<CouplingManager>
  build(std::string_view type, std::string_view scope);

  COUPLING_TYPE type_ = COUPLING_TYPE::NONE;
  COUPLING_SCOPE scope_ = COUPLING_SCOPE::NONE;
  std::array<CouplingInterface, 2u> interfaces_;
  std::array<std::unique_ptr<MeshCoupling>, 2u> meshes_;
  std::array<std::unordered_map<std::string, std::unique_ptr<FieldCoupling>>, 2u>
      fields_;
};

} // namespace cocoa
