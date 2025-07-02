#pragma once

#include "plugins.hpp"

#ifdef COCOA_ENABLE_MEDCOUPLING

// medcoupling
#include <MEDCouplingRemapper.hxx>

// local
#include "coupling/coupling_manager.hpp"
#include "coupling/field_med.hpp"

namespace cocoa
{

struct Problem;

struct CouplingMED: public CouplingManager
{
  CouplingMED(COUPLING_SCOPE scope): CouplingManager(COUPLING_TYPE::MEDCOUPLING, scope)
  {}
  ~CouplingMED() = default;

  void setup(CouplingInterface interfaceSrc, CouplingInterface interfaceTgt) override;
  void initFieldCoupling() override;

  void project(std::string_view fieldSrc, std::string_view fieldTgt) override;

  MEDCoupling::MEDCouplingRemapper remapper_;
  FieldMED mask_;
  // INTERP_KERNEL::IntersectionType interpType_ = INTERP_KERNEL::Geometric2D;
  INTERP_KERNEL::IntersectionType interpType_ = INTERP_KERNEL::Triangulation;
};

} // namespace cocoa

#endif
