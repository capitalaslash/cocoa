// std
#include <cassert>

// local
#include "coupling/coupling_manager.hpp"
#include "enums.hpp"
#include "problem/problem.hpp"

int main()
{
  using namespace cocoa;

  auto p1 = Problem::build(PROBLEM_TYPE::FD1D, EQN_TYPE::NONE);
  auto p2 = Problem::build("fd1d", "none");

  p1->setup({{"config_file", std::filesystem::path{"fd1d_heat.dat"}}});
  p2->setup({{"config_file", std::filesystem::path{"fd1d_hc.dat"}}});
  p1->couplingType_ = COUPLING_TYPE::SIMPLE;
  p2->couplingType_ = COUPLING_TYPE::SIMPLE;

  // check that both problems have been set to use the same coupling type
  assert(p1->couplingType_ == p2->couplingType_);

  auto coupling12 = CouplingManager::build(p1->couplingType_);
  coupling12->setup(p1.get(), p2.get());

  // MEDManager coupling21;
  // coupling21.setup(p2.get(), p1.get());

  p1->print();
  p2->print();

  while (p1->run() || p2->run())
  {
    // coupling21.project("T", "Tsys");
    p1->advance();
    p1->solve();
    p1->print();

    coupling12->project("T", "uExternal");
    p2->advance();
    p2->solve();
    p2->print();
  }

  return 0;
}
