#pragma once

// std
#include <cstdlib>
#include <string_view>

// coupling enum =======================================================
enum struct COUPLING_TYPE : int8_t
{
  NONE = 0,
  SIMPLE = 1,
  MEDCOUPLING = 2,
  OFM2M = 3,
};

inline COUPLING_TYPE str2coupling(std::string_view name)
{
  if (name == "simple")
    return COUPLING_TYPE::SIMPLE;
  if (name == "medcoupling")
    return COUPLING_TYPE::MEDCOUPLING;
  std::abort();
  return COUPLING_TYPE::NONE;
}

// problem enum ========================================================
enum struct PROBLEM_TYPE : int8_t
{
  NONE = 0,
  FD1D = 1,
  FEMUS = 2,
  PROXPDE = 3,
  OFCOM = 4,
  OFORG = 5,
};

inline PROBLEM_TYPE str2problem(std::string_view name)
{
  if (name == "fd1d")
    return PROBLEM_TYPE::FD1D;
  if (name == "femus")
    return PROBLEM_TYPE::FEMUS;
  if (name == "proxpde")
    return PROBLEM_TYPE::PROXPDE;
  if (name == "oforg")
    return PROBLEM_TYPE::OFORG;
  std::abort();
  return PROBLEM_TYPE::NONE;
}

// fdbc enum ===========================================================
enum struct FDBC_TYPE : int8_t
{
  NONE = 0,
  DIRICHLET = 1,
  NEUMANN = 2,
};

inline FDBC_TYPE str2fdbc(std::string_view buffer)
{
  if (buffer == "dirichlet")
    return FDBC_TYPE::DIRICHLET;
  if (buffer == "neumann")
    return FDBC_TYPE::NEUMANN;
  std::abort();
  return FDBC_TYPE::NONE;
}
