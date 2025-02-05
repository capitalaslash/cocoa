#include "problem/fdutils.hpp"

// =====================================================================
template <>
bool ParamsFD::set<uint>(std::string_view name, uint value)
{
  auto const & [_, success] = data_.emplace(
      name,
      std::pair{FD_PARAM_TYPE::INTEGER, std::tuple{value, 0.0, std::vector<double>{}}});
  return success;
}

template <>
bool ParamsFD::set<double>(std::string_view name, double value)
{
  auto const & [_, success] = data_.emplace(
      name,
      std::pair{FD_PARAM_TYPE::SCALAR, std::tuple{0u, value, std::vector<double>{}}});
  return success;
}

template <>
bool ParamsFD::set<std::vector<double>>(
    std::string_view name, std::vector<double> value)
{
  auto const & [_, success] =
      data_.emplace(name, std::pair{FD_PARAM_TYPE::VECTOR, std::tuple{0u, 0.0, value}});
  return success;
}
