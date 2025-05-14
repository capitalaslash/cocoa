#pragma once

// std
#include <array>
#include <cstdint>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <vector>

// libfmt
#include <fmt/core.h>
#include <fmt/ranges.h>

// local
#include "la.hpp"

namespace cocoa
{

// =====================================================================
enum struct FD_BC_TYPE : uint8_t
{
  NONE = 0,
  DIRICHLET,
  NEUMANN,
};

inline FD_BC_TYPE str2FDBCType(std::string_view name)
{
  if (name == "dirichlet")
    return FD_BC_TYPE::DIRICHLET;
  if (name == "neumann")
    return FD_BC_TYPE::NEUMANN;
  fmt::print(stderr, "boundary condition type not recognized: {}\n", name);
  std::abort();
  return FD_BC_TYPE::NONE;
}

enum struct FD_BC_SIDE : uint8_t
{
  LEFT = 0u,
  RIGHT = 1u,
  BOTTOM = 2u,
  TOP = 3u,
  FRONT = 4u,
  BACK = 5u,
};

// =====================================================================
struct FDBC
{
  FD_BC_SIDE side;
  FD_BC_TYPE type = FD_BC_TYPE::NONE;
  VectorFD values;
  VectorFD ghostValues;

  FDBC() = default;
  ~FDBC() = default;

  FDBC(FD_BC_SIDE s, FD_BC_TYPE const t, double const value, size_t const size = 1u):
      side{s},
      type{t},
      values(size, value),
      ghostValues(size, 0.0)
  {}

  FDBC(FD_BC_SIDE s, FD_BC_TYPE const t, VectorFD const & v):
      side{s},
      type{t},
      values{v},
      ghostValues(v.size(), 0.0)
  {}

  auto
  init(FD_BC_SIDE s, FD_BC_TYPE const t, double const value, size_t const size = 1u)
      -> void
  {
    side = s;
    type = t;
    values.data_.resize(size, value);
    ghostValues.data_.resize(size, 0.0);
  }

  auto init(FD_BC_SIDE s, FD_BC_TYPE const t, VectorFD const & v) -> void
  {
    side = s;
    type = t;
    values.data_ = v.data_;
    ghostValues.data_.resize(v.size(), 0.0);
  }
};

template <uint8_t dim>
struct FDBCList
{
  std::array<FDBC, 2 * dim> data_;
  auto left() -> FDBC & { return data_[0]; }
  auto right() -> FDBC & { return data_[1]; }

  auto bottom() -> FDBC &
  {
    static_assert(dim >= 2u, "bottom() is not available in 1d");
    return data_[2];
  }
  auto top() -> FDBC &
  {
    static_assert(dim >= 2u, "top() is not available in 1d");
    return data_[3];
  }
};

using FDBCList1D = FDBCList<1u>;
using FDBCList2D = FDBCList<2u>;

// =====================================================================
template <uint8_t dim>
struct MeshFD
{
  using Int_T = std::array<uint, dim>;
  using Real_T = std::array<double, dim>;

  Real_T start_;
  Int_T n_;
  Real_T h_;

  MeshFD() = default;
  MeshFD(Real_T const & start, Real_T const & end, Int_T const & n)
  {
    init(start, end, n);
  }
  ~MeshFD() = default;

  auto init(Real_T const & start, Real_T end, Int_T n)
  {
    start_ = start;
    for (uint d = 0u; d < dim; d++)
    {
      n_[d] = n[d] + 1;
      h_[d] = (end[d] - start[d]) / n[d];
    }
  }

  constexpr auto nPts() const -> uint
  {
    auto n = 1u;
    for (uint d = 0u; d < dim; d++)
      n *= n_[d];
    return n;
  }

  constexpr auto nElems() const -> uint
  {
    auto n = 1u;
    for (uint d = 0u; d < dim; d++)
      n *= n_[d] - 1;
    return n;
  }

  constexpr auto end() const -> Real_T
  {
    Real_T end;
    for (uint d = 0u; d < dim; d++)
      end[d] = start_[d] + (n_[d] - 1) * h_[d];
    return end;
  }

  auto pt(Int_T const & ijk) const -> Real_T
  {
    Real_T p;
    for (uint d = 0u; d < dim; d++)
    {
      assert(ijk[d] < n_[d]);
      p[d] = start_[d] + ijk[d] * h_[d];
    }
    return p;
  }

  // auto pt(uint const i) const -> double { return pt({i})[0]; }
};

using MeshFD1D = MeshFD<1u>;
using MeshFD2D = MeshFD<2u>;

// =====================================================================
enum class FD_PARAM_TYPE : std::uint8_t
{
  INTEGER = 0u,
  SCALAR = 1u,
  VECTOR = 2u,
};

inline FD_PARAM_TYPE str2FDParamType(std::string_view name)
{
  using enum FD_PARAM_TYPE;
  if (name == "integer")
    return INTEGER;
  else if (name == "scalar")
    return SCALAR;
  else if (name == "vector")
    return VECTOR;
  fmt::print(stderr, "param type {} not recognized\n", name);
  std::abort();
  return INTEGER;
}

inline std::string fdParamType2str(FD_PARAM_TYPE const t)
{
  switch (t)
  {
    using enum FD_PARAM_TYPE;
  case INTEGER:
    return "integer";
  case SCALAR:
    return "scalar";
  case VECTOR:
    return "vector";
  default:
    fmt::print(
        stderr,
        "FDParamType {} not recognized\n",
        static_cast<std::underlying_type_t<FD_PARAM_TYPE>>(t));
    std::abort();
  }
}

struct ParamsFD
{
  using Param_T = std::tuple<uint, double, std::vector<double>>;
  using Value_T = std::pair<FD_PARAM_TYPE, Param_T>;
  using ParamList_T = std::unordered_map<std::string, Value_T>;

  template <typename T>
  auto set(std::string_view name, T value) -> bool;

  template <FD_PARAM_TYPE T>
  auto get(std::string_view name) const
  {
    // TODO: assert that the type corresponds
    if (!data_.contains(std::string(name)))
    {
      fmt::print(stderr, "parameter {} not available\n", name);
      std::abort();
    }
    using UT = std::underlying_type_t<FD_PARAM_TYPE>;
    return std::get<static_cast<UT>(T)>(data_.at(std::string(name)).second);
  }

  ParamList_T data_;
};

} // namespace cocoa

// formatters ==========================================================

template <>
struct fmt::formatter<cocoa::FDBC>
{
  constexpr auto parse(format_parse_context & ctx) -> format_parse_context::iterator
  {
    return ctx.begin();
  }

  auto format(cocoa::FDBC const bc, format_context & ctx) const
      -> format_context::iterator
  {
    return fmt::format_to(
        ctx.out(),
        "s: {} t: {} v: {} g: {}",
        static_cast<uint8_t>(bc.side),
        static_cast<uint8_t>(bc.type),
        bc.values,
        bc.ghostValues);
  }
};

template <>
struct fmt::formatter<cocoa::FD_PARAM_TYPE>
{
  constexpr auto parse(format_parse_context & ctx) -> format_parse_context::iterator
  {
    return ctx.begin();
  }

  auto format(cocoa::FD_PARAM_TYPE const t, format_context & ctx) const
      -> format_context::iterator
  {
    return fmt::format_to(ctx.out(), "{}", fdParamType2str(t));
  }
};

template <>
struct fmt::formatter<cocoa::ParamsFD>
{
  constexpr auto parse(format_parse_context & ctx) -> format_parse_context::iterator
  {
    return ctx.begin();
  }

  auto format(cocoa::ParamsFD const & params, format_context & ctx) const
      -> format_context::iterator
  {
    std::string data = "\n";
    for (auto const & [name, value]: params.data_)
    {
      data += fmt::format("  {}: ", name);
      switch (value.first)
      {
        using enum cocoa::FD_PARAM_TYPE;
      case INTEGER:
      {
        auto const value = params.get<INTEGER>(name);
        data += fmt::format("{}\n", value);
        break;
      }
      case SCALAR:
      {
        auto const value = params.get<SCALAR>(name);
        data += fmt::format("{:e}\n", value);
        break;
      }
      case VECTOR:
      {
        auto const value = params.get<VECTOR>(name);
        data += fmt::format("{::e}\n", value);
        break;
      }
      }
    }
    return fmt::format_to(ctx.out(), "{}", data);
  }
};
