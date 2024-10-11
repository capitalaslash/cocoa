#pragma once

// std
#include <array>
#include <cstdint>
#include <string_view>
#include <tuple>
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
struct FDBC
{
  FD_BC_TYPE type = FD_BC_TYPE::NONE;
  double value = 0.0;
};

// =====================================================================
enum struct FD_SOLVER_TYPE : uint8_t
{
  NONE = 0,
  TRIDIAG,
  VANKA1D,
  VANKA2DCB,
  VANKA2DSCI,
};

inline FD_SOLVER_TYPE str2fdsolver(std::string_view name)
{
  if (name == "tridiag")
    return FD_SOLVER_TYPE::TRIDIAG;
  if (name == "vanka1d")
    return FD_SOLVER_TYPE::VANKA1D;
  if (name == "vanka2dcb")
    return FD_SOLVER_TYPE::VANKA2DCB;
  if (name == "vanka2dsci")
    return FD_SOLVER_TYPE::VANKA2DSCI;
  fmt::print(stderr, "solver {} not recognized\n", name);
  std::abort();
  return FD_SOLVER_TYPE::NONE;
}

// =====================================================================
// struct VectorFD
// {
//   VectorFD() = default;
//   explicit VectorFD(size_t const n): data_(n) {}
//   ~VectorFD() = default;

//   double * data() { return data_.data(); }

//   std::vector<double> data_;
// };
using VectorFD = std::vector<double>;

// =====================================================================
struct MatrixTriDiag
{
  MatrixTriDiag() = default;
  explicit MatrixTriDiag(size_t const n): diag(n), diagUp(n), diagDown(n) {}
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
struct MatrixCSR
{
  struct Entry
  {
    size_t clm;
    double value;
  };

  using Triplet_T = std::tuple<size_t, size_t, double>;
  using Row_T = std::vector<Entry>;

  MatrixCSR() = default;
  explicit MatrixCSR(size_t n): n_{n}, data_{n_, Row_T{}} {}
  ~MatrixCSR() = default;

  void init(size_t n);

  void clear();
  void close();

  size_t n_;
  std::vector<Triplet_T> triplets_;
  std::vector<Row_T> data_;
};

std::vector<double> operator*(MatrixCSR const & m, std::vector<double> const & v);

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
  using Vec2D_T = std::array<double, 2U>;
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

  Vec2D_T operator*(Vec2D_T const & v) const
  {
    Vec2D_T result;
    result[0] = data[0][0] * v[0] + data[0][1] * v[1];
    result[1] = data[1][0] * v[0] + data[1][1] * v[1];
    return result;
  }

  std::array<Vec2D_T, 2U> data;
};
