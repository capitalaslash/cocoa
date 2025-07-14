#pragma once

// std
#include <filesystem>
#include <memory>
#include <string_view>
#include <vector>

// fmt
#include <fmt/core.h>

// local
#include "enums.hpp"

namespace cocoa
{

struct MeshCoupling
{
  MeshCoupling() = default;
  explicit MeshCoupling(COUPLING_TYPE type, COUPLING_SCOPE scope)
      : type_{type}
      , scope_{scope}
  {}
  virtual ~MeshCoupling() = default;

  virtual void init(
      std::string_view name,
      COUPLING_SCOPE const scope,
      Marker const marker,
      std::string_view bdName,
      uint const dim,
      uint const spaceDim,
      std::vector<double> const & coords,
      std::vector<uint> const & conn,
      std::vector<uint> const & offsets) = 0;

  virtual void printVTK(std::filesystem::path const & path) = 0;

  static std::unique_ptr<MeshCoupling> build(COUPLING_TYPE type);
  static std::unique_ptr<MeshCoupling> build(std::string_view type);

  COUPLING_TYPE type_ = COUPLING_TYPE::NONE;
  COUPLING_SCOPE scope_ = COUPLING_SCOPE::NONE;
  Marker marker_;
  std::string bdName_ = "";
};

struct MeshSimple: public MeshCoupling
{
  MeshSimple(): MeshCoupling{COUPLING_TYPE::SIMPLE, COUPLING_SCOPE::NONE} {}
  ~MeshSimple() = default;

  void init(
      std::string_view name,
      COUPLING_SCOPE const scope,
      Marker marker,
      std::string_view bdName,
      uint const dim,
      uint const spaceDim,
      std::vector<double> const & coords,
      std::vector<uint> const & conn,
      std::vector<uint> const & offsets) override
  {
    scope_ = scope;
    marker_ = marker;
    bdName_ = bdName;
  }

  void printVTK(std::filesystem::path const & /*path*/) override
  {
    if (messageVTK_)
    {
      fmt::println(stderr, "VTK print is not implemented in MeshSimple.");
      messageVTK_ = false;
    }
  }

  bool messageVTK_ = true;
};

} // namespace cocoa
