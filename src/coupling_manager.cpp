#include "coupling_manager.hpp"

void CouplingManager::setup(Problem * pSrc, Problem * pTgt)
{
  pSrc_ = pSrc;
  pTgt_ = pTgt;
}

void CouplingManager::project(
    std::string_view fieldNameSrc, std::string_view fieldNameTgt)
{}