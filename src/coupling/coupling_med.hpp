#pragma once

#include "plugins.hpp"

#ifdef COCOA_ENABLE_MEDCOUPLING

// medcoupling
#include <MEDCouplingRemapper.hxx>

// local
#include "coupling/coupling_manager.hpp"

namespace cocoa
{

struct Problem;

struct CouplingMED: public CouplingManager
{
  CouplingMED(): CouplingManager(COUPLING_TYPE::MEDCOUPLING) {}
  ~CouplingMED() = default;

  void setup(Problem * src, Problem * tgt);
  void project(std::string_view srcName, std::string_view tgtName);

  Problem * problemSrc_;
  Problem * problemTgt_;
  MEDCoupling::MEDCouplingRemapper remapper;
  // INTERP_KERNEL::IntersectionType interpType_ = INTERP_KERNEL::Geometric2D;
  INTERP_KERNEL::IntersectionType interpType_ = INTERP_KERNEL::Triangulation;
};

} // namespace cocoa

#endif
