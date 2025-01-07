#pragma once

// std
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
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

inline double norm2sq(VectorFD const & v)
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

  std::vector<std::pair<uint, double>> row(uint k) const
  {
    // left
    if (k == 0u)
      return {{k, diags_[1][k]}, {k + 1, diags_[2][k]}};
    // right
    else if (k == diags_[0].size() - 1)
      return {{k - 1, diags_[0][k]}, {k, diags_[1][k]}};
    // inside
    return {{k - 1, diags_[0][k]}, {k, diags_[1][k]}, {k + 1, diags_[2][k]}};
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

  // std::vector<double> & diag() { return diags_[1]; }
  // std::vector<double> & diagUp() { return diags_[2]; }
  // std::vector<double> & diagDown() { return diags_[0]; }
  // std::vector<double> const & diag() const { return diags_[1]; }
  // std::vector<double> const & diagUp() const { return diags_[2]; }
  // std::vector<double> const & diagDown() const { return diags_[0]; }

  std::array<std::vector<double>, 3U> diags_;
};

std::vector<double> operator*(MatrixTriDiag const & m, std::vector<double> const & v);

template <>
struct fmt::formatter<MatrixTriDiag>
{
  template <typename ParseContext>
  constexpr auto parse(ParseContext & ctx);

  template <typename FormatContext>
  auto format(MatrixTriDiag const & m, FormatContext & ctx);
};

template <typename ParseContext>
constexpr auto fmt::formatter<MatrixTriDiag>::parse(ParseContext & ctx)
{
  return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<MatrixTriDiag>::format(MatrixTriDiag const & m, FormatContext & ctx)
{
  return fmt::format_to(
      ctx.out(), "down: {}\ndiag: {}\nup: {}", m.diags_[0], m.diags_[1], m.diags_[2]);
}

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

  Row_T const & row(uint k) const { return data_[k]; }

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

template <>
struct fmt::formatter<MatrixCSR::Entry>
{
  template <typename ParseContext>
  constexpr auto parse(ParseContext & ctx) -> format_parse_context::iterator;

  template <typename FormatContext>
  auto format(MatrixCSR::Entry const & e, FormatContext & ctx) const
      -> format_context::iterator;
};

template <typename ParseContext>
constexpr auto fmt::formatter<MatrixCSR::Entry>::parse(ParseContext & ctx)
    -> format_parse_context::iterator
{
  return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<MatrixCSR::Entry>::format(
    MatrixCSR::Entry const & e, FormatContext & ctx) const -> format_context::iterator
{
  return fmt::format_to(ctx.out(), "[{}, {}]", e.clm, e.value);
}

template <>
struct fmt::formatter<MatrixCSR>
{
  template <typename ParseContext>
  constexpr auto parse(ParseContext & ctx);

  template <typename FormatContext>
  auto format(MatrixCSR const & m, FormatContext & ctx);
};

template <typename ParseContext>
constexpr auto fmt::formatter<MatrixCSR>::parse(ParseContext & ctx)
{
  return ctx.begin();
}

template <typename FormatContext>
auto fmt::formatter<MatrixCSR>::format(MatrixCSR const & m, FormatContext & ctx)
{
  return fmt::format_to(
      ctx.out(), "stored entries:\n{}\npending entries:\n{}", m.data_, m.triplets_);
}

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

// =====================================================================
template <typename Matrix>
double computeResidual(
    Matrix const & m,
    std::vector<double> const & x,
    std::vector<double> const & b,
    double const area)
{
  std::vector<double> res(x.size());
  auto const tmp = m * x;
  for (uint k = 0; k < res.size(); k++)
    res[k] = (b[k] - tmp[k]);
  double const resNorm = std::sqrt(norm2sq(res) * area);
  return resNorm;
}

struct SolverInfo
{
  uint nIters;
  double residual;
};

// template <typename Matrix>
// struct Solver
// {
//   Solver() = default;
//   virtual ~Solver() = default;

//   auto set(Matrix const & m) -> void { m_ = &m; }
//   auto virtual solve(VectorFD const & b, VectorFD & x) const -> SolverInfo = 0;

//   std::unique_ptr<Matrix const> m_;
//   uint maxIters = 1000u;
//   double tolerance = 1e-6;
// };

// template <typename Matrix>
// using Solver_T =
//     std::function<SolverInfo(Solver<Matrix> const *, VectorFD const &, VectorFD
//     &)>;

// struct SolverTriDiag: public Solver<MatrixTriDiag>
// {
//   auto virtual solve(VectorFD const & b, VectorFD & x) const
//       -> SolverInfo override final;
// };

template <typename Matrix>
using Solver_T = std::function<SolverInfo(
    Matrix const &, VectorFD const &, VectorFD &, double const, uint const)>;

SolverInfo solveTriDiag(MatrixTriDiag const & m, VectorFD const & b, VectorFD & x);

SolverInfo solveVanka1D(
    MatrixTriDiag const & m,
    VectorFD const & b,
    VectorFD & x,
    double const tolerance,
    uint maxIters);

template <typename Matrix>
SolverInfo solveJacobi(
    Matrix const & m,
    VectorFD const & b,
    VectorFD & x,
    double const tolerance,
    uint maxIters)
{
  uint const n = b.size();

  std::vector<double> uNew(x.size());
  double const rhsNorm = std::sqrt(norm2sq(b) / n);
  fmt::print("rhsNorm: {:.8e}\n", rhsNorm);

  for (uint i = 0; i < maxIters; i++)
  {
    double const resNorm = computeResidual(m, x, b, 1.0 / n);
    // fmt::print("iter: {:3d}, current residual: {:.8e}\n", j, resNorm);
    if (resNorm < tolerance * rhsNorm)
    {
      return {i, resNorm / rhsNorm};
    }

    for (uint k = 0U; k < x.size(); k++)
    {
      double valueNew = b[k];
      double diag = 0.0;
      for (auto const & [clm, value]: m.row(k))
      {
        if (clm == k)
        {
          diag = value;
        }
        else
        {
          valueNew -= value * x[clm];
        }
      }
      uNew[k] = valueNew / diag;
    }

    // update current solution
    for (uint k = 0U; k < uNew.size(); k++)
      x[k] = uNew[k];
  }

  return {maxIters, computeResidual(m, x, b, 1.0 / n) / rhsNorm};
}

template <typename Matrix>
SolverInfo solveGaussSeidel(
    Matrix const & m,
    VectorFD const & b,
    VectorFD & x,
    double const tolerance,
    uint maxIters)
{
  uint const n = b.size();

  double const rhsNorm = std::sqrt(norm2sq(b) / n);
  fmt::print("rhsNorm: {:.8e}\n", rhsNorm);

  for (uint i = 0; i < maxIters; i++)
  {
    double const resNorm = computeResidual(m, x, b, 1.0 / n);
    // fmt::print("iter: {:3d}, current residual: {:.8e}\n", j, resNorm);
    if (resNorm < tolerance * rhsNorm)
    {
      return {i, resNorm / rhsNorm};
    }

    for (uint k = 0U; k < x.size(); k++)
    {
      double valueNew = b[k];
      double diag = 0.0;
      for (auto const [clm, value]: m.data_[k])
      {
        if (clm == k)
        {
          diag = value;
        }
        else
        {
          valueNew -= value * x[clm];
        }
      }
      x[k] = valueNew / diag;
    }
  }

  return {maxIters, computeResidual(m, x, b, 1.0 / n) / rhsNorm};
}

inline double solveLine(
    MatrixCSR const & m,
    std::vector<double> const & rhs,
    std::vector<double> const & u,
    size_t const id)
{
  assert(m.data_[id][0].clm == id);
  double value = rhs[id];
  for (uint j = 1; j < m.data_[id].size(); j++)
  {
    value -= m.data_[id][j].value * u[m.data_[id][j].clm];
  }
  return value / m.data_[id][0].value;
}

template <typename Matrix>
SolverInfo solveVanka2DCB(
    Matrix const & m,
    VectorFD const & b,
    VectorFD & x,
    double const tolerance,
    uint maxIters)
{
  uint const n = b.size();

  double const rhsNorm = std::sqrt(norm2sq(b) / n);
  fmt::print("rhsNorm: {:.8e}\n", rhsNorm);

  for (uint i = 0U; i < maxIters; i++)
  {
    double const resNorm = computeResidual(m, x, b, 1.0 / n);

    // fmt::print("iter: {:3d}, current residual: {:.8e}\n", j, resNorm);
    if (resNorm < tolerance * rhsNorm)
    {
      return {i, resNorm / rhsNorm};
    }

    // checkerboard
    for (uint k = 0U; k < n; k += 2U)
      x[k] = solveLine(m, b, x, k);
    for (uint k = 1U; k < n; k += 2U)
      x[k] = solveLine(m, b, x, k);
  }

  return {maxIters, computeResidual(m, x, b, 1.0 / n) / rhsNorm};
}

template <typename Matrix>
SolverInfo solveVanka2DSCI(
    Matrix const & /*m*/,
    VectorFD const & /*b*/,
    VectorFD & /*x*/,
    double const /*tolerance*/,
    uint /*maxIters*/)
{
  return {0u, 0.0};

  // uint const n = b.size();

  // double const rhsNorm = std::sqrt(norm2sq(b)  / n);
  // fmt::print("rhsNorm: {:.8e}\n", rhsNorm);

  // for (uint i = 0U; i < maxIters; i++)
  // {
  //   double const resNorm = computeResidual(m, x, b, 1.0 / n);

  //   // fmt::print("iter: {:3d}, current residual: {:.8e}\n", j, resNorm);
  //   if (resNorm < tolerance * rhsNorm)
  //   {
  //     return {i, resNorm / rhsNorm};
  //   }

  //   // for (uint k = 0U; k < x.size(); k++)
  //   //   uOld_[k] = x[k];

  //   // sides + corners + inside
  //   // sides
  //   for (uint k = 0U; k < 4U; k++)
  //   {
  //     auto const dofList = sideDOF(n, k);
  //     for (auto const & dof: dofList)
  //       x[dof] = solveLine(m, b, x, dof);
  //   }

  //   // corners
  //   for (auto const & id: cornerDOF(n))
  //     x[id] = solveLine(m, b, x, id);

  //   // inside
  //   for (uint j = 1U; j < n[1] - 1; j += 1U)
  //     for (uint i = 1U; i < n[0] - 1; i += 1U)
  //     {
  //       auto const id = i + j * n[0];
  //       x[id] = solveLine(m, b, x, id);
  //     }
  // }

  // return {maxIters, computeResidual(m, x, b, 1.0 / n) / rhsNorm};
}
