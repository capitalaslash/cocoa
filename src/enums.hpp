#pragma once

// std
#include <cstdlib>
#include <string>
#include <string_view>

// local
#include "plugins.hpp"

// coupling enum =======================================================
enum struct COUPLING_TYPE : int8_t
{
  NONE = 0,
  SIMPLE,
  MEDCOUPLING,
  OFM2M,
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
  FD1D,
  FD2D,
  FEMUS,
  PROXPDE,
  OFCOM,
  OFORG,
};

inline PROBLEM_TYPE str2problem(std::string_view name)
{
  if (name == "fd1d")
    return PROBLEM_TYPE::FD1D;
  if (name == "fd2d")
    return PROBLEM_TYPE::FD2D;
  if (name == "femus")
    return PROBLEM_TYPE::FEMUS;
  if (name == "proxpde")
    return PROBLEM_TYPE::PROXPDE;
  if (name == "oforg")
    return PROBLEM_TYPE::OFORG;
  std::abort();
  return PROBLEM_TYPE::NONE;
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
  CUSTOM,
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

// MED_CELL_TYPE =======================================================
#ifdef COCOA_ENABLE_MEDCOUPLING

#include <MCIdType.hxx>
#include <NormalizedGeometricTypes> // missing extension!!!

enum struct MED_CELL_TYPE : int8_t
{
  LINE2 = INTERP_KERNEL::NORM_SEG2,
  TRIANGLE3 = INTERP_KERNEL::NORM_TRI3,
  QUAD4 = INTERP_KERNEL::NORM_QUAD4,
  HEX8 = INTERP_KERNEL::NORM_HEXA8,
};

inline mcIdType MEDCellTypeToIKCell(MED_CELL_TYPE t)
{
  return static_cast<mcIdType>(t);
}

#else

using mcIdType = std::uint32_t;

enum struct MED_CELL_TYPE : int8_t
{
  LINE2 = 12,
  TRIANGLE3 = 23,
  QUAD4 = 24,
  HEX8 = 38,
};

inline mcIdType MEDCellTypeToIKCell(MED_CELL_TYPE t)
{
  return static_cast<mcIdType>(t);
}

#endif
