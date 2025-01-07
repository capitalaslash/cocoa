#pragma once

// std
#include <algorithm>
#include <array>
#include <cassert>
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
  TRIDIAG,
  VANKA1D,
  JACOBI2D,
  GAUSSSEIDEL2D,
  VANKA2DCB,
  VANKA2DSCI,
};

inline FD_SOLVER_TYPE str2fdsolver(std::string_view name)
{
  if (name == "tridiag")
    return FD_SOLVER_TYPE::TRIDIAG;
  if (name == "vanka1d")
    return FD_SOLVER_TYPE::VANKA1D;
  if (name == "jacobi2d")
    return FD_SOLVER_TYPE::JACOBI2D;
  if (name == "gaussseidel2d")
    return FD_SOLVER_TYPE::GAUSSSEIDEL2D;
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

inline double norm2sq(std::vector<double> & v)
{
  double sum = 0.0;
  for (auto const & value: v)
  {
    sum += value * value;
  }
  return sum;
}

// =====================================================================
struct MatrixTriDiag
{
  MatrixTriDiag() = default;
  explicit MatrixTriDiag(size_t const n):
      diags_{std::vector<double>(n), std::vector<double>(n), std::vector<double>(n)}
  {}
  ~MatrixTriDiag() = default;

  void init(size_t n)
  {
    diags_[0].resize(n, 0.0);
    diags_[1].resize(n, 0.0);
    diags_[2].resize(n, 0.0);
  }

  auto size() const -> size_t { return diags_[1].size(); }

  void set(int const row, int const clm, double const value)
  {
    assert(std::abs(clm - row) <= 1);
    diags_[clm - row + 1][row] = value;
  }

  void add(int const row, int const clm, double const value)
  {
    assert(std::abs(clm - row) <= 1);
    diags_[clm - row + 1][row] += value;
  }

  auto at(uint const row, uint const clm) const -> double
  {
    if (clm == row - 1)
      return diags_[0][row];
    else if (clm == row)
      return diags_[1][row];
    else if (clm == row + 1)
      return diags_[2][row];
    else
      return 0.0;
  }

  void clearRow(uint const row)
  {
    diags_[0][row] = 0.0;
    diags_[1][row] = 0.0;
    diags_[2][row] = 0.0;
  }
  void clear()
  {
    uint const n = diags_[1].size();
    std::vector<double>(n).swap(diags_[0]);
    std::vector<double>(n).swap(diags_[1]);
    std::vector<double>(n).swap(diags_[2]);
  }
  void close() {}

  std::vector<double> & diag() { return diags_[1]; }
  std::vector<double> & diagUp() { return diags_[2]; }
  std::vector<double> & diagDown() { return diags_[0]; }
  std::vector<double> const & diag() const { return diags_[1]; }
  std::vector<double> const & diagUp() const { return diags_[2]; }
  std::vector<double> const & diagDown() const { return diags_[0]; }

  std::array<std::vector<double>, 3U> diags_;
};

std::vector<double> operator*(MatrixTriDiag const & m, std::vector<double> const & v);

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

  auto size() const -> size_t { return n_; }

  void set(int const /*row*/, int const /*clm*/, double const /*value*/)
  {
    fmt::print(stderr, "set() is not supported in MatrixCSR!\n");
    std::abort();
  }

  void add(uint const row, uint const clm, double const value)
  {
    triplets_.emplace_back(row, clm, value);
  }

  auto at(uint const row, uint const clm) const -> double
  {
    for (auto const & entry: data_[row])
      if (entry.clm == clm)
        return entry.value;
    return 0.0;
  }

  void clearRow(uint const row);
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
