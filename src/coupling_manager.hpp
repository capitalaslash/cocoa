#pragma once

// std
#include <string_view>

struct Problem;

struct CouplingManager
{
  CouplingManager() = default;
  virtual ~CouplingManager() = default;

  virtual void setup(Problem * pSrc, Problem * pTgt) = 0;

  virtual void
  project(std::string_view fieldNameSrc, std::string_view fieldNameTgt) = 0;

  Problem * pSrc_;
  Problem * pTgt_;
};
