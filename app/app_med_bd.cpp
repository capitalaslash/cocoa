// local
#include "coupling/coupling_manager.hpp"
#include "problem/problem.hpp"

int main()
{
  using namespace cocoa;

  auto p1 = Problem::build(PROBLEM_TYPE::PROXPDE, EQN_TYPE::NS);
  auto p2 = Problem::build(PROBLEM_TYPE::PROXPDE, EQN_TYPE::NS);

  p1->setup({{"config_file", std::filesystem::path{"proxpde_med_bd1.yaml"}}});
  p2->setup({{"config_file", std::filesystem::path{"proxpde_med_bd2.yaml"}}});

  p1->print();
  p2->print();

  auto coupling12 =
      CouplingManager::build(COUPLING_TYPE::MEDCOUPLING, COUPLING_SCOPE::BOUNDARY);
  coupling12->setup(
      {p1.get(), "top", {"vel"}},
      {p2.get(), "bottom", {"vel"}},
      INTERPOLATION_METHOD::P1P1);
  auto coupling21 =
      CouplingManager::build(COUPLING_TYPE::MEDCOUPLING, COUPLING_SCOPE::BOUNDARY);
  coupling21->setup(
      {p2.get(), "bottom", {"p"}},
      {p1.get(), "top", {"p"}},
      INTERPOLATION_METHOD::P1P1);

  while (p1->run() || p2->run())
  {
    p1->advance();
    p2->advance();

    for (auto k = 0u; k < 2u; k++)
    {
      coupling21->project("p");
      p1->solve();
      coupling12->project("vel");
      p2->solve();
    }

    p1->print();
    p2->print();
  }

  return 0;
}
