// stl
#include <memory>

// cocoa
#include <cocoa/coupling/field_coupling.hpp>
#include <cocoa/coupling/field_med.hpp>
#include <cocoa/enums.hpp>
#include <cocoa/plugins.hpp>

namespace cocoa
{

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
#else
  case COUPLING_TYPE::MEDCOUPLING:
  {
    fmt::println(stderr, "MED coupling not available, reverting to Simple coupling");
    return std::unique_ptr<FieldCoupling>{new FieldSimple};
    break;
  }
#endif
#ifdef COCOA_ENABLE_OFM2M
  case COUPLING_TYPE::OFM2M:
  {
    return std::unique_ptr<FieldCoupling>{new FieldOForg};
    break;
  }
#else
  case COUPLING_TYPE::OFM2M:
  {
    fmt::println(stderr, "OFM2M coupling not available, reverting to Simple coupling");
    return std::unique_ptr<FieldCoupling>{new FieldSimple};
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
  return FieldCoupling::build(str2couplingType(couplingType));
}

} // namespace cocoa
