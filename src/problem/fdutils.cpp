#include "problem/fdutils.hpp"

// std
#include <cassert>

// =====================================================================
std::vector<double> operator*(MatrixTriDiag const & m, std::vector<double> const & v)
{
  uint const n = m.diag().size();
  assert(n == v.size());
  std::vector<double> r(n);

  r[0] = m.diag()[0] * v[0] + m.diagUp()[0] * v[1];
  for (uint k = 1; k < n - 1; k++)
  {
    r[k] = m.diagDown()[k] * v[k - 1] + m.diag()[k] * v[k] + m.diagUp()[k] * v[k + 1];
  }
  r[n - 1] = m.diagDown()[n - 1] * v[n - 2] + m.diag()[n - 1] * v[n - 1];

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

void MatrixCSR::clear()
{
  std::vector<Triplet_T>{}.swap(triplets_);
  std::vector<Row_T>(n_).swap(data_);
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
