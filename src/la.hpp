#pragma once

// stl
#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <functional>
#include <vector>

// fmtlib
#include <fmt/core.h>
#include <fmt/ranges.h>

// =====================================================================
struct VectorFD
{
  VectorFD() = default;
  ~VectorFD() = default;
  explicit VectorFD(size_t const n, double const value = 0.0): data_(n, value) {}
  auto resize(size_t const n, double const value = 0.0) -> void
  {
    data_.resize(n, value);
  }
  explicit VectorFD(std::vector<double> const & data): data_{data} {}

  auto size() const -> size_t { return data_.size(); }
  auto data() -> double * { return data_.data(); }
  auto operator[](uint const k) const -> double { return data_[k]; }

  auto add(uint const k, double const value) -> void
  {
    assert(k < data_.size());
    data_[k] += value;
  }
  auto set(uint const k, double const value) -> void
  {
    assert(k < data_.size());
    data_[k] = value;
  }
  auto setRange(uint const start, uint const end, double const value)
  {
    assert(end - start <= data_.size());
    std::fill(data_.begin() + start, data_.begin() + end, value);
  }
  auto zero() -> void { std::fill(data_.begin(), data_.end(), 0.0); }
  auto norm2sq() const -> double;

  std::vector<double> data_;
};

template <>
struct fmt::formatter<VectorFD>
{
  constexpr auto parse(format_parse_context & ctx) -> format_parse_context::iterator
  {
    return ctx.begin();
  }

  auto format(VectorFD const & v, format_context & ctx) const
      -> format_context::iterator
  {
    return fmt::format_to(ctx.out(), "{::e}", v.data_);
  }
};

// =====================================================================
struct MatrixTriDiag
{
  MatrixTriDiag() = default;
  explicit MatrixTriDiag(size_t const n): diags_{VectorFD{n}, VectorFD{n}, VectorFD{n}}
  {}
  ~MatrixTriDiag() = default;

  void init(size_t n)
  {
    diags_[0].data_.resize(n, 0.0);
    diags_[1].data_.resize(n, 0.0);
    diags_[2].data_.resize(n, 0.0);
  }

  auto size() const -> size_t { return diags_[1].size(); }

  void set(int const row, int const clm, double const value)
  {
    assert(std::abs(clm - row) <= 1);
    diags_[clm - row + 1].add(row, value);
  }

  void add(int const row, int const clm, double const value)
  {
    assert(std::abs(clm - row) <= 1);
    diags_[clm - row + 1].add(row, value);
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
    diags_[0].set(row, 0.0);
    diags_[1].set(row, 0.0);
    diags_[2].set(row, 0.0);
  }
  void clear()
  {
    uint const n = diags_[1].size();
    std::vector<double>(n).swap(diags_[0].data_);
    std::vector<double>(n).swap(diags_[1].data_);
    std::vector<double>(n).swap(diags_[2].data_);
  }
  void close() {}

  std::array<VectorFD, 3U> diags_;
};

auto operator*(MatrixTriDiag const & m, VectorFD const & v) -> VectorFD;

template <>
struct fmt::formatter<MatrixTriDiag>
{
  constexpr auto parse(format_parse_context & ctx) -> format_parse_context::iterator
  {
    return ctx.begin();
  }

  auto format(MatrixTriDiag const & m, format_context & ctx) const
      -> format_context::iterator
  {
    return fmt::format_to(
        ctx.out(), "down: {}\ndiag: {}\nup: {}", m.diags_[0], m.diags_[1], m.diags_[2]);
  }
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
  MatrixCSR(size_t n, size_t nnz);
  ~MatrixCSR() = default;

  void init(size_t n, size_t nnz);

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

  void print_sparsity_pattern(std::filesystem::path const & path);

  size_t n_;
  size_t nnz_ = 0u;
  bool isClosed_ = false;
  std::vector<Triplet_T> triplets_;
  std::vector<Row_T> data_;
};

auto operator*(MatrixCSR const & m, VectorFD const & v) -> VectorFD;

template <>
struct fmt::formatter<MatrixCSR::Entry>
{
  constexpr auto parse(format_parse_context & ctx) -> format_parse_context::iterator
  {
    return ctx.begin();
  }

  auto format(MatrixCSR::Entry const & e, format_context & ctx) const
      -> format_context::iterator
  {
    return fmt::format_to(ctx.out(), "[{}, {}]", e.clm, e.value);
  }
};

template <>
struct fmt::formatter<MatrixCSR>
{
  constexpr auto parse(format_parse_context & ctx) -> format_parse_context::iterator
  {
    return ctx.begin();
  }

  auto format(MatrixCSR const & m, format_context & ctx) const
      -> format_context::iterator
  {
    return fmt::format_to(
        ctx.out(), "stored entries:\n{}\npending entries:\n{}", m.data_, m.triplets_);
  }
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
    Matrix const & m, VectorFD const & x, VectorFD const & b, double const area)
{
  VectorFD res(x.size());
  auto const tmp = m * x;
  for (uint k = 0; k < res.size(); k++)
    res.set(k, b[k] - tmp[k]);
  double const resNorm = std::sqrt(res.norm2sq() * area);
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

  VectorFD xNew(x.size());
  double const rhsNorm = std::sqrt(b.norm2sq() / n);
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
      xNew.set(k, valueNew / diag);
    }

    // update current solution
    x = xNew;
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

  double const rhsNorm = std::sqrt(b.norm2sq() / n);
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
      x.set(k, valueNew / diag);
    }
  }

  return {maxIters, computeResidual(m, x, b, 1.0 / n) / rhsNorm};
}

inline double solveLine(
    MatrixCSR const & m, VectorFD const & rhs, VectorFD const & u, size_t const id)
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

  double const rhsNorm = std::sqrt(b.norm2sq() / n);
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
      x.set(k, solveLine(m, b, x, k));
    for (uint k = 1U; k < n; k += 2U)
      x.set(k, solveLine(m, b, x, k));
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
  fmt::print(stderr, "solveVanka2DSCI not implemented\n");
  std::abort();
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
