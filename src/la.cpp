#include "la.hpp"

// =====================================================================
std::vector<double> operator*(MatrixTriDiag const & m, std::vector<double> const & v)
{
  uint const n = m.size();
  assert(n == v.size());
  std::vector<double> r(n);

  r[0] = m.at(0, 0) * v[0] + m.at(0, 1) * v[1];
  for (uint k = 1u; k < n - 1; k++)
  {
    r[k] = m.at(k, k - 1) * v[k - 1] + m.at(k, k) * v[k] + m.at(k, k + 1) * v[k + 1];
  }
  r[n - 1] = m.at(n - 1, n - 2) * v[n - 2] + m.at(n - 1, n - 1) * v[n - 1];

  return r;
}

// =====================================================================
void MatrixCSR::init(size_t n)
{
  n_ = n;
  data_.resize(n, Row_T{});
  for (uint k = 0; k < n; k++)
    data_[k].reserve(5);
}

void MatrixCSR::clearRow(uint const row)
{
  // clear store data
  Row_T{}.swap(data_[row]);
  data_[row].reserve(5);

  // clear triplet data
  triplets_.erase(
      std::remove_if(
          triplets_.begin(),
          triplets_.end(),
          [row](auto const & triplet)
          {
            if (std::get<0>(triplet) == row)
              return true;
            return false;
          }),
      triplets_.end());
}

void MatrixCSR::clear()
{
  std::vector<Triplet_T>{}.swap(triplets_);
  std::vector<Row_T>(n_).swap(data_);
  isClosed_ = false;
}

void MatrixCSR::close()
{
  for (auto const & [row, clm, value]: triplets_)
  {
    assert(row < n_ && clm < n_);
    auto clmFound = false;
    for (uint k = 0; k < data_[row].size(); k++)
    {
      // check if clm has already been set
      if (data_[row][k].clm == clm)
      {
        clmFound = true;
        data_[row][k].value += value;
        break;
      }
    }
    // when the clm is not yet stored, create new storage place
    if (!clmFound)
    {
      data_[row].emplace_back(clm, value);
    }
  }

  // clear out triplets after storing them
  std::vector<Triplet_T>().swap(triplets_);
  isClosed_ = true;
}

void MatrixCSR::print_sparsity_pattern(std::filesystem::path const & path)
{
  assert(isClosed_);
  std::FILE * f = std::fopen(path.string().c_str(), "w");

  for (auto const & row: data_)
  {
    std::string line(n_ + 2, ' ');
    line[0] = '|';
    line[line.size() - 1] = '|';
    for (auto const [clm, _]: row)
      line[clm + 1] = 'o';
    fmt::print(f, "{}\n", line);
  }

  std::fclose(f);
}

std::vector<double> operator*(MatrixCSR const & m, std::vector<double> const & v)
{
  assert(m.n_ == v.size());
  std::vector<double> r(m.n_, 0.0);
  for (uint k = 0U; k < m.n_; k++)
  {
    auto const & row = m.data_[k];
    for (auto const & [id, value]: row)
    {
      r[k] += value * v[id];
    }
  }

  return r;
}

// =====================================================================
SolverInfo solveTriDiag(MatrixTriDiag const & m, VectorFD const & b, VectorFD & x)
{
  uint const n = b.size();
  double const rhsNorm = std::sqrt(norm2sq(b) / n);
  fmt::print("rhsNorm: {:.8e}\n", rhsNorm);

  std::vector<double> upPrime(n);
  std::vector<double> rhsPrime(n);

  upPrime[0] = m.at(0, 1) / m.at(0, 0);
  for (uint k = 1U; k < n - 1; k++)
  {
    upPrime[k] = m.at(k, k + 1) / (m.at(k, k) - m.at(k, k - 1) * upPrime[k - 1]);
  }

  rhsPrime[0] = b[0] / m.at(0, 0);
  for (uint k = 1U; k < n; k++)
  {
    rhsPrime[k] = (b[k] - m.at(k, k - 1) * rhsPrime[k - 1]) /
                  (m.at(k, k) - m.at(k, k - 1) * upPrime[k - 1]);
  }

  x[n - 1] = rhsPrime[n - 1];
  for (int k = n - 2; k >= 0; k--)
  {
    x[k] = rhsPrime[k] - upPrime[k] * x[k + 1];
  }

  return {1U, computeResidual(m, x, b, 1.0 / n) / rhsNorm};
}

SolverInfo solveVanka1D(
    MatrixTriDiag const & m,
    VectorFD const & b,
    VectorFD & x,
    double const tolerance,
    uint const maxIters)
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

    for (uint k = 1U; k < n - 1; k += 2U)
    {
      x[k] =
          (b[k] - m.at(k, k - 1) * x[k - 1] - m.at(k, k + 1) * x[k + 1]) / m.at(k, k);
    }
    for (uint k = 2U; k < n - 1; k += 2U)
    {
      x[k] =
          (b[k] - m.at(k, k - 1) * x[k - 1] - m.at(k, k + 1) * x[k + 1]) / m.at(k, k);
    }
    x[0] = (b[0] - m.at(0, 1) * x[1]) / m.at(0, 0);
    x[n - 1] = (b[n - 1] - m.at(n - 1, n - 2) * x[n - 2]) / m.at(n - 1, n - 1);
  }

  return {maxIters, computeResidual(m, x, b, 1.0 / n) / rhsNorm};
}
