// std
#include <memory>

// local
#include "plugins.hpp"

#include "enums.hpp"
#include "field_coupling.hpp"
#include "field_med.hpp"

std::unique_ptr<FieldCoupling> FieldCoupling::build(COUPLING_TYPE type)
{
  switch (type)
  {
  case COUPLING_TYPE::SIMPLE:
  {
    return std::unique_ptr<FieldCoupling>{new FieldSimple};
    break;
  }
#ifdef COCOA_ENABLE_MEDCOUPLING
  case COUPLING_TYPE::MEDCOUPLING:
  {
    return std::unique_ptr<FieldCoupling>{new FieldMED};
    break;
  }
#endif
#ifdef COCOA_ENABLE_OFM2M
  case COUPLING_TYPE::OFM2M:
  {
    return std::unique_ptr<FieldCoupling>{new FieldOForg};
    break;
  }
#endif
  default:
  {
    std::abort();
  }
  }
}

std::unique_ptr<FieldCoupling> FieldCoupling::build(std::string_view couplingType)
{
  return FieldCoupling::build(str2coupling(couplingType));
}
