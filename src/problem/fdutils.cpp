#include "problem/fdutils.hpp"

// std
#include <cassert>

void MatrixCSR::init(size_t n)
{
  n_ = n;
  data_.resize(n, Row_T{});
  for (uint k = 0; k < n; k++)
    data_[k].reserve(5);
}

void MatrixCSR::clear()
{
  std::vector<Triplet_T>().swap(triplets_);
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
      if (data_[row][k].clm == clm)
      {
        clmFound = true;
        data_[row][k].value += value;
        break;
      }
    }
    if (!clmFound)
    {
      data_[row].emplace_back(clm, value);
    }
  }
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
