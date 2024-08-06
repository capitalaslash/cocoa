#pragma once

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

inline Problem * buildProblem(PROBLEM_TYPE type)
{
  switch (type)
  {
  case PROBLEM_TYPE::FD1D:
  {
    return new ProblemFD1D{};
    break;
  }
  case PROBLEM_TYPE::PROXPDE:
  {
    return new ProblemProXPDE{};
    break;
  }
  default:
  {
    std::abort();
  }
  }
}
