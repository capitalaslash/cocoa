// std
#include <memory>

// fmtlib
#include <fmt/core.h>

// local
#include "plugins.hpp"

#include "coupling/coupling_manager.hpp"
#include "coupling/coupling_med.hpp"
#include "coupling/coupling_simple.hpp"
#include "enums.hpp"

namespace cocoa
{

std::unique_ptr<CouplingManager> CouplingManager::build(COUPLING_TYPE type)
{
  switch (type)
  {
  case COUPLING_TYPE::SIMPLE:
  {
    return std::unique_ptr<CouplingManager>{new CouplingSimple};
    break;
  }
#ifdef COCOA_ENABLE_MEDCOUPLING
  case COUPLING_TYPE::MEDCOUPLING:
  {
    return std::unique_ptr<CouplingManager>{new CouplingMED};
    break;
  }
#else
  case COUPLING_TYPE::MEDCOUPLING:
  {
    fmt::println(stderr, "MED coupling not available, reverting to Simple coupling");
    return std::unique_ptr<CouplingManager>{new CouplingSimple};
    break;
  }
#endif
#ifdef COCOA_ENABLE_OFM2M
  case COUPLING_TYPE::OFM2M:
  {
    return std::unique_ptr<CouplingManager>{new CouplingOFM2M};
    break;
  }
#else
  case COUPLING_TYPE::OFM2M:
  {
    fmt::println(stderr, "OFM2M coupling not available, reverting to Simple coupling");
    return std::unique_ptr<CouplingManager>{new CouplingSimple};
    break;
  }
#endif
  default:
  {
    fmt::println(stderr, "coupling type not recognized!");
    std::abort();
  }
  }
}

inline std::unique_ptr<CouplingManager>
CouplingManager::build(std::string_view couplingType)
{
  return CouplingManager::build(str2coupling(couplingType));
}

} // namespace cocoa
