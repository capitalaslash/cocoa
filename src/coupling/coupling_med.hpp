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

  void setup(
      CouplingInterface interfaceSrc,
      CouplingInterface interfaceTgt,
      INTERPOLATION_METHOD method) override;
  void initFieldCoupling() override;

  void project(std::string_view fieldSrc, std::string_view fieldTgt) override;

  MEDCoupling::MEDCouplingRemapper remapper_;
  FieldMED mask_;
  // optons: Geometric2D, Triangulation, Convex, PointLocator
  INTERP_KERNEL::IntersectionType interpolationType_ = INTERP_KERNEL::Triangulation;
};

} // namespace cocoa

#endif
