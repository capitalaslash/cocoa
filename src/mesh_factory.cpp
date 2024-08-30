// std
#include <memory>

// local
#include "plugins.hpp"

#include "enums.hpp"
#include "mesh_coupling.hpp"
#include "mesh_med.hpp"

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
#endif
#ifdef COCOA_ENABLE_OFM2M
  case COUPLING_TYPE::OFM2M:
  {
    return std::unique_ptr<MeshCoupling>{new MeshOForg};
    break;
  }
#endif
  default:
  {
    std::abort();
  }
  }
}

std::unique_ptr<MeshCoupling> MeshCoupling::build(std::string_view couplingType)
{
  return MeshCoupling::build(str2coupling(couplingType));
}
