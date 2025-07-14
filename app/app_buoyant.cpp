// local
#include "coupling/coupling_manager.hpp"
#include "enums.hpp"
#include "problem/problem.hpp"

int main()
{
  using namespace cocoa;

  auto pNS = Problem::build(PROBLEM_TYPE::PROXPDE, EQN_TYPE::NS_BOUSSINESQ);
  pNS->setup({{"config_file", std::filesystem::path{"proxpde_buoyant_ns.yaml"}}});

  auto pHeat = Problem::build(PROBLEM_TYPE::PROXPDE, EQN_TYPE::HEAT_BUOYANT);
  pHeat->setup({{"config_file", std::filesystem::path{"proxpde_buoyant_heat.yaml"}}});

  auto couplingHeatToNS =
      CouplingManager::build(COUPLING_TYPE::MEDCOUPLING, COUPLING_SCOPE::VOLUME);
  couplingHeatToNS->setup(
      {pHeat.get(), {"T"}}, {pNS.get(), {"T"}}, INTERPOLATION_METHOD::P1P1);

  auto couplingNSToHeat =
      CouplingManager::build(COUPLING_TYPE::MEDCOUPLING, COUPLING_SCOPE::VOLUME);
  couplingNSToHeat->setup(
      {pNS.get(), {"vel"}}, {pHeat.get(), {"vel"}}, INTERPOLATION_METHOD::P1P1);

  pNS->print();
  pHeat->print();

  while (pNS->run() || pHeat->run())
  {
    couplingHeatToNS->project("T", "T");
    pNS->advance();
    pNS->solve();
    pNS->print();

    couplingNSToHeat->project("vel", "vel");
    pHeat->advance();
    pHeat->solve();
    pHeat->print();
  }

  return 0;
}
