#pragma once

#include "plugins.hpp"

// stl
#include <cstdlib>
#include <string>
#include <string_view>

// fmt
#include "fmt/core.h"

// medcoupling
#ifdef COCOA_ENABLE_MEDCOUPLING
#include <MCIdType.hxx>
#include <NormalizedGeometricTypes> // missing extension!!!
#endif

namespace cocoa
{

// coupling enum =======================================================
enum struct COUPLING_TYPE : uint8_t
{
  NONE = 0u,
  SIMPLE,
  MEDCOUPLING,
  OFM2M,
};

inline COUPLING_TYPE str2coupling(std::string_view name)
{
  if (name == "simple")
    return COUPLING_TYPE::SIMPLE;
  else if (name == "medcoupling")
    return COUPLING_TYPE::MEDCOUPLING;
  else if (name == "ofm2m")
    return COUPLING_TYPE::OFM2M;
  else
  {
    fmt::println(stderr, "Coupling type {} not recognized!", name);
    std::abort();
  }
  return COUPLING_TYPE::NONE;
}

inline std::string couplingType2str(COUPLING_TYPE const type)
{
  using enum COUPLING_TYPE;
  switch (type)
  {
  case SIMPLE:
    return "simple";
  case MEDCOUPLING:
    return "medcoupling";
  case OFM2M:
    return "ofm2m";
  default:
    return "none";
  }
}

// problem enum ========================================================
enum struct PROBLEM_TYPE : uint8_t
{
  NONE = 0u,
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
enum struct EQN_TYPE : uint8_t
{
  NONE = 0u,
  HEAT,
  HEAT_COUPLED,
  HEAT_BUOYANT,
  HEAT_OC,
  NS,
  NS_BOUSSINESQ,
  CUSTOM,
};

inline EQN_TYPE str2eqn(std::string_view name)
{
  using enum EQN_TYPE;
  if (name == "heat")
    return HEAT;
  else if (name == "heat_coupled")
    return HEAT_COUPLED;
  else if (name == "heat_buoyant")
    return HEAT_BUOYANT;
  else if (name == "heat_oc")
    return HEAT_OC;
  else if (name == "ns")
    return NS;
  else if (name == "ns_boussinesq")
    return NS_BOUSSINESQ;
  else if (name == "none")
    return NONE;
  fmt::println(stderr, "equation type {} not recognized", name);
  std::abort();
  return NONE;
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
    return "heat_coupled";
  case EQN_TYPE::HEAT_BUOYANT:
    return "heat_buoyant";
  case EQN_TYPE::HEAT_OC:
    return "heat_oc";
  case EQN_TYPE::NS:
    return "ns";
  case EQN_TYPE::NS_BOUSSINESQ:
    return "ns_boussinesq";
  default:
    fmt::print(
        stderr,
        "equation type {} not recognized\n",
        static_cast<std::underlying_type_t<EQN_TYPE>>(type));
    std::abort();
    return "";
  }
}

// MED_CELL_TYPE =======================================================
#ifdef COCOA_ENABLE_MEDCOUPLING

enum struct MED_CELL_TYPE : uint8_t
{
  LINE2 = INTERP_KERNEL::NormalizedCellType::NORM_SEG2,
  TRIANGLE3 = INTERP_KERNEL::NormalizedCellType::NORM_TRI3,
  QUAD4 = INTERP_KERNEL::NormalizedCellType::NORM_QUAD4,
  HEX8 = INTERP_KERNEL::NormalizedCellType::NORM_HEXA8,
};

inline mcIdType MEDCellTypeToIKCell(MED_CELL_TYPE t)
{
  return static_cast<mcIdType>(t);
}

#else

using mcIdType = std::uint32_t;

enum struct MED_CELL_TYPE : uint8_t
{
  LINE2 = 12u,
  TRIANGLE3 = 23u,
  QUAD4 = 24u,
  HEX8 = 38u,
};

inline mcIdType MEDCellTypeToIKCell(MED_CELL_TYPE t)
{
  return static_cast<mcIdType>(t);
}

#endif

} // namespace cocoa