#pragma once

// std
#include <string_view>
#include <vector>

// local
#include "enums.hpp"

struct MeshCoupling
{
  MeshCoupling() = default;
  explicit MeshCoupling(COUPLING_TYPE type): type_{type} {}
  virtual ~MeshCoupling() = default;

  virtual void init(
      std::string_view name,
      uint const dim,
      uint const spaceDim,
      std::vector<double> const & coords,
      std::vector<uint> conn,
      std::vector<uint> offsets) = 0;

  COUPLING_TYPE type_ = COUPLING_TYPE::NONE;
};

struct MeshSimple: public MeshCoupling
{
  MeshSimple(): MeshCoupling{COUPLING_TYPE::SIMPLE} {}
  ~MeshSimple() = default;

  void init(
      std::string_view name,
      uint const dim,
      uint const spaceDim,
      std::vector<double> const & coords,
      std::vector<uint> conn,
      std::vector<uint> offsets) override
  {}
};
