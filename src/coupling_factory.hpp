#pragma once

#include "defs.hpp"

// std
#include <memory>

// local
#include "coupling_manager.hpp"
#include "coupling_simple.hpp"

enum struct COUPLING_TYPE : char
{
  NONE = 0,
  SIMPLE = 1,
  MED = 2,
  OFM2M = 3,
};

inline std::unique_ptr<CouplingManager> buildCoupling(COUPLING_TYPE type)
{
  switch (type)
  {
  case COUPLING_TYPE::SIMPLE:
  {
    return std::unique_ptr<CouplingManager>{new CouplingSimple{}};
    break;
  }
#ifdef COCOA_ENABLE_MED
  case COUPLING_TYPE::MED:
  {
    return std::unique_ptr<CouplingManager>{new CouplingMED{}};
    break;
  }
#endif
#ifdef COCOA_ENABLE_OFM2M
  case COUPLING_TYPE::OFM2M:
  {
    return std::unique_ptr<CouplingManager>{new CouplingOFM2M{}};
    break;
  }
#endif
  default:
  {
    std::abort();
  }
  }
}