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
  explicit MeshCoupling(COUPLING_TYPE type): type_{type} {}
  virtual ~MeshCoupling() = default;

  virtual void init(
      std::string_view name,
      uint const dim,
      uint const spaceDim,
      std::vector<double> const & coords,
      std::vector<uint> conn,
      std::vector<uint> offsets) = 0;

  virtual void printVTK(std::filesystem::path const & path) = 0;

  static std::unique_ptr<MeshCoupling> build(COUPLING_TYPE type);
  static std::unique_ptr<MeshCoupling> build(std::string_view type);

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
