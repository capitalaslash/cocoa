// std
#include <cassert>

// local
#include "coupling/coupling_manager.hpp"
#include "problem/problem.hpp"
#include "problem/problem_proxpde.hpp"

int main()
{
  using namespace cocoa;

  auto p1 = Problem::build(PROBLEM_TYPE::PROXPDE, EQN_TYPE::HEAT);
  p1->setup({{"config_file", std::filesystem::path{"proxpde_heat.yaml"}}});
  dynamic_cast<ProblemProXPDEHeat *>(p1.get())->fieldsP0_["q"]
      << [](proxpde::Vec3 const & p) { return std::sin(M_PI * p[1]); };

  auto p2 = Problem::build("proxpde", "heat_coupled");
  p2->setup({{"config_file", std::filesystem::path{"proxpde_hc.yaml"}}});

  // check that both problems have been set to use the same coupling type
  assert(p1->couplingType_ == p2->couplingType_);

  auto coupling12 = CouplingManager::build(p1->couplingType_);
  coupling12->setup(p1.get(), p2.get());

  p1->print();
  p2->print();

  while (p1->run() || p2->run())
  {
    p1->advance();
    p1->solve();
    p1->print();

    coupling12->project("Tcfd");
    p2->advance();
    p2->solve();
    p2->print();
  }

  return 0;
}
