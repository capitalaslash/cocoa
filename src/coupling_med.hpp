#pragma once

// medcoupling
#include <MEDCouplingRemapper.hxx>

// local
#include "coupling_manager.hpp"

struct Problem;

struct CouplingMED: public CouplingManager
{
  CouplingMED() = default;
  ~CouplingMED() = default;

  void setup(Problem * src, Problem * tgt);
  void project(std::string_view srcName, std::string_view tgtName);

  Problem * problemSrc_;
  Problem * problemTgt_;
  MEDCoupling::MEDCouplingRemapper remapper;
};
