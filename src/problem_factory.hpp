#pragma once

// std
#include <memory>

// local
#include "plugins.hpp"

#include "problem.hpp"
#include "problem_fd1d.hpp"
#include "problem_femus.hpp"
#include "problem_oforg.hpp"
#include "problem_proxpde.hpp"

enum struct PROBLEM_TYPE : char
{
  NONE = 0,
  FD1D = 1,
  FEMUS = 2,
  PROXPDE = 3,
  OFCOM = 4,
  OFORG = 5,
};

inline std::unique_ptr<Problem> buildProblem(PROBLEM_TYPE type)
{
  switch (type)
  {
  case PROBLEM_TYPE::FD1D:
  {
    return std::unique_ptr<Problem>{new ProblemFD1D};
    break;
  }
#ifdef COCOA_ENABLE_FEMUS
  case PROBLEM_TYPE::FEMUS:
  {
    return std::unique_ptr<Problem>{new ProblemFEMUS};
    break;
  }
#endif
#ifdef COCOA_ENABLE_OFORG
  case PROBLEM_TYPE::OFORG:
  {
    return std::unique_ptr<Problem>{new ProblemOForg};
    break;
  }
#endif
#ifdef COCOA_ENABLE_PROXPDE
  case PROBLEM_TYPE::PROXPDE:
  {
    return std::unique_ptr<Problem>{new ProblemProXPDE};
    break;
  }
#endif
  default:
  {
    std::abort();
  }
  }
}

