// std
#include <memory>

// fmt
#include <fmt/core.h>

// local
#include "plugins.hpp"

#include "coupling/mesh_coupling.hpp"
#include "coupling/mesh_med.hpp"
#include "enums.hpp"

namespace cocoa
{

std::unique_ptr<MeshCoupling> MeshCoupling::build(COUPLING_TYPE type)
{
  switch (type)
  {
  case COUPLING_TYPE::SIMPLE:
  {
    return std::unique_ptr<MeshCoupling>{new MeshSimple};
    break;
  }
#ifdef COCOA_ENABLE_MEDCOUPLING
  case COUPLING_TYPE::MEDCOUPLING:
  {
    return std::unique_ptr<MeshCoupling>{new MeshMED};
    break;
  }
#else
  case COUPLING_TYPE::MEDCOUPLING:
  {
    fmt::print(stderr, "MED coupling not available, reverting to Simple coupling\n");
    return std::unique_ptr<MeshCoupling>{new MeshSimple};
    break;
  }
#endif
#ifdef COCOA_ENABLE_OFM2M
  case COUPLING_TYPE::OFM2M:
  {
    return std::unique_ptr<MeshCoupling>{new MeshOForg};
    break;
  }
#else
  case COUPLING_TYPE::OFM2M:
  {
    fmt::print(stderr, "OFM2M coupling not available, reverting to Simple coupling\n");
    return std::unique_ptr<MeshCoupling>{new MeshSimple};
    break;
  }
#endif
  default:
  {
    fmt::print(stderr, "coupling type has not been set.\n");
    std::abort();
  }
  }
}

std::unique_ptr<MeshCoupling> MeshCoupling::build(std::string_view couplingType)
{
  return MeshCoupling::build(str2coupling(couplingType));
}

} // namespace cocoa
