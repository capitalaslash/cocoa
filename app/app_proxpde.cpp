// std
#include <cassert>

// local
#include "coupling/coupling_manager.hpp"
#include "enums.hpp"
#include "problem/problem.hpp"
#include "problem/problem_proxpde.hpp"

int main()
{
  using namespace cocoa;

  auto p1 = Problem::build(PROBLEM_TYPE::PROXPDE, EQN_TYPE::HEAT);
  p1->setup({{"config_file", std::filesystem::path{"proxpde_heat.yaml"}}});
  dynamic_cast<ProblemProXPDEHeat *>(p1.get())->fieldsP0_["q"]
      << [](proxpde::Vec3 const & p) { return std::sin(M_PI * p[1]); };
  p1->print();

  auto p2 = Problem::build("proxpde", "heat_coupled");
  p2->setup({{"config_file", std::filesystem::path{"proxpde_hc.yaml"}}});
  p2->print();

  auto coupling12 =
      CouplingManager::build(COUPLING_TYPE::MEDCOUPLING, COUPLING_SCOPE::VOLUME);
  coupling12->setup({p1.get(), {"T"}}, {p2.get(), {"Tcfd"}}, INTERPOLATION_METHOD::P1P1);

  while (p1->run() || p2->run())
  {
    p1->advance();
    p1->solve();
    p1->print();

    coupling12->project("T", "Tcfd");
    p2->advance();
    p2->solve();
    p2->print();
  }

  return 0;
}
