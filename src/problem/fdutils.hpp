#pragma once

// std
#include <array>
#include <cstdint>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <vector>

// fmt
#include <fmt/core.h>
#include <fmt/ranges.h>

// =====================================================================
enum struct FD_BC_TYPE : int8_t
{
  NONE = 0,
  DIRICHLET,
  NEUMANN,
};

inline FD_BC_TYPE str2fdbc(std::string_view name)
{
  if (name == "dirichlet")
    return FD_BC_TYPE::DIRICHLET;
  if (name == "neumann")
    return FD_BC_TYPE::NEUMANN;
  fmt::print(stderr, "boundary condition type not recognized: {}\n", name);
  std::abort();
  return FD_BC_TYPE::NONE;
}

// =====================================================================
struct FDBC
{
  FD_BC_TYPE type = FD_BC_TYPE::NONE;
  std::vector<double> values;
};

// =====================================================================
enum struct FD_SOLVER_TYPE : uint8_t
{
  NONE = 0,
  GAUSSSEIDEL,
  JACOBI,
  TRIDIAG,
  VANKA1D,
  VANKA2DCB,
  VANKA2DSCI,
};

inline FD_SOLVER_TYPE str2fdsolver(std::string_view name)
{
  if (name == "gaussseidel")
    return FD_SOLVER_TYPE::GAUSSSEIDEL;
  else if (name == "jacobi")
    return FD_SOLVER_TYPE::JACOBI;
  else if (name == "tridiag")
    return FD_SOLVER_TYPE::TRIDIAG;
  else if (name == "vanka1d")
    return FD_SOLVER_TYPE::VANKA1D;
  else if (name == "vanka2dcb")
    return FD_SOLVER_TYPE::VANKA2DCB;
  else if (name == "vanka2dsci")
    return FD_SOLVER_TYPE::VANKA2DSCI;
  else
  {
    fmt::print(stderr, "solver {} not recognized\n", name);
    std::abort();
  }
  return FD_SOLVER_TYPE::NONE;
}
