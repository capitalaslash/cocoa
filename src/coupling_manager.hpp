#pragma once

// std
#include <string_view>

struct Problem;

struct CouplingManager
{
  CouplingManager() = default;
  virtual ~CouplingManager() = default;

  virtual void setup(Problem * pSrc, Problem * pTgt);

  virtual void project(std::string_view fieldNameSrc, std::string_view fieldNameTgt);

  Problem * pSrc_;
  Problem * pTgt_;
};
