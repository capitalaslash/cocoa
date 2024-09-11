// std
#include <memory>

// local
#include "plugins.hpp"

#include "problem.hpp"
#include "problem_fd1d.hpp"
#include "problem_femus.hpp"
#include "problem_oforg.hpp"
#include "problem_proxpde.hpp"

std::unique_ptr<Problem> Problem::build(PROBLEM_TYPE problemType, EQN_TYPE eqnType)
{
  switch (problemType)
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
    return std::unique_ptr<Problem>{ProblemProXPDE::build(eqnType)};
    break;
  }
#endif
  default:
  {
    std::abort();
  }
  }
}

std::unique_ptr<Problem>
Problem::build(std::string_view problemType, std::string_view eqnType)
{
  return Problem::build(str2problem(problemType), str2eqn(eqnType));
}