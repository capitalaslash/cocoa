#include "problem/fdutils.hpp"

// std
#include <cassert>

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
