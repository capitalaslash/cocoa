#pragma once

// std
#include <array>
#include <cstdint>
#include <string_view>
#include <vector>

// fmt
#include <fmt/core.h>

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
  std::abort();
  return FD_BC_TYPE::NONE;
}

// =====================================================================
enum struct FD_SOLVER_TYPE : uint8_t
{
  NONE = 0,
  TRIDIAG,
  VANKA,
};

inline FD_SOLVER_TYPE str2fdsolver(std::string_view name)
{
  if (name == "tridiag")
    return FD_SOLVER_TYPE::TRIDIAG;
  if (name == "vanka")
    return FD_SOLVER_TYPE::VANKA;
  fmt::print(stderr, "solver {} not recognized\n", name);
  std::abort();
  return FD_SOLVER_TYPE::NONE;
}

// =====================================================================
struct MatrixTriDiag
{
  MatrixTriDiag() = default;
  explicit MatrixTriDiag(size_t n): diag(n), diagUp(n), diagDown(n) {}
  ~MatrixTriDiag() = default;

  void init(size_t n)
  {
    diag.resize(n);
    diagUp.resize(n);
    diagDown.resize(n);
  }

  std::vector<double> diag;
  std::vector<double> diagUp;
  std::vector<double> diagDown;
};

// =====================================================================
template <size_t n>
struct MatrixDense
{
  MatrixDense() = default;
  ~MatrixDense() = default;

  double & operator()(size_t i, size_t j) { return data[i][j]; }
  double operator()(size_t i, size_t j) const { return data[i][j]; }

  std::array<std::array<double, n>, n> data;
};

template <>
struct MatrixDense<2U>
{
  MatrixDense() = default;
  ~MatrixDense() = default;

  double & operator()(size_t i, size_t j) { return data[i][j]; }
  double operator()(size_t i, size_t j) const { return data[i][j]; }

  double det() const { return data[0][0] * data[1][1] - data[1][0] * data[0][1]; }
  MatrixDense<2U> inverse() const
  {
    double const idet = 1. / det();
    MatrixDense<2U> im;
    im(0, 0) = data[1][1] * idet;
    im(0, 1) = -data[0][1] * idet;
    im(1, 0) = -data[1][0] * idet;
    im(1, 1) = data[0][0] * idet;
    return im;
  }

  std::array<double, 2U> operator*(std::array<double, 2U> const & v) const
  {
    std::array<double, 2U> result;
    result[0] = data[0][0] * v[0] + data[0][1] * v[1];
    result[1] = data[1][0] * v[0] + data[1][1] * v[1];
    return result;
  }

  std::array<std::array<double, 2U>, 2U> data;
};
