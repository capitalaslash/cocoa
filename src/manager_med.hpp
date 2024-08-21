#pragma once

// medcoupling
#include <MEDCouplingRemapper.hxx>

struct Problem;

struct MEDManager
{
  MEDManager() = default;
  ~MEDManager() = default;

  void setup(Problem * src, Problem * tgt);
  void project(std::string_view srcName, std::string_view tgtName);

  Problem * problemSrc_;
  Problem * problemTgt_;
  MEDCoupling::MEDCouplingRemapper remapper;
};
