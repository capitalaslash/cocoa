#pragma once

// std
#include <cstdlib>
#include <string>
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

inline FDBC_TYPE str2fdbc(std::string_view name)
{
  if (name == "dirichlet")
    return FDBC_TYPE::DIRICHLET;
  if (name == "neumann")
    return FDBC_TYPE::NEUMANN;
  std::abort();
  return FDBC_TYPE::NONE;
}

// equation type enum ===============================================
enum struct EQN_TYPE : int8_t
{
  NONE = 0,
  HEAT,
  HEAT_COUPLED,
  HEAT_BUOYANT,
  NS,
  NS_BOUSSINESQ,
};

inline EQN_TYPE str2eqn(std::string_view name)
{
  if (name == "heat")
    return EQN_TYPE::HEAT;
  else if (name == "heatCoupled")
    return EQN_TYPE::HEAT_COUPLED;
  else if (name == "heatBuoyant")
    return EQN_TYPE::HEAT_BUOYANT;
  else if (name == "ns")
    return EQN_TYPE::NS;
  else if (name == "nsBoussinesq")
    return EQN_TYPE::NS_BOUSSINESQ;
  else if (name == "none")
    return EQN_TYPE::NONE;
  abort();
  return EQN_TYPE::NONE;
}

inline std::string eqn2str(EQN_TYPE type)
{
  switch (type)
  {
  case EQN_TYPE::NONE:
    return "none";
  case EQN_TYPE::HEAT:
    return "heat";
  case EQN_TYPE::HEAT_COUPLED:
    return "heatCoupled";
  case EQN_TYPE::HEAT_BUOYANT:
    return "heatBuoyant";
  case EQN_TYPE::NS:
    return "ns";
  case EQN_TYPE::NS_BOUSSINESQ:
    return "nsBoussinesq";
  default:
    abort();
    return "";
  }
}
