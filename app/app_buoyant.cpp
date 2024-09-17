// local
#include "coupling/coupling_manager.hpp"
#include "problem/problem.hpp"

int main()
{
  auto pNS = Problem::build(PROBLEM_TYPE::PROXPDE, EQN_TYPE::NS_BOUSSINESQ);
  pNS->setup({{"config_file", "proxpde_buoyant_ns.yaml"}});

  auto pHeat = Problem::build(PROBLEM_TYPE::PROXPDE, EQN_TYPE::HEAT_BUOYANT);
  pHeat->setup({{"config_file", "proxpde_buoyant_heat.yaml"}});

  auto couplingHeatToNS = CouplingManager::build(COUPLING_TYPE::MEDCOUPLING);
  couplingHeatToNS->setup(pHeat.get(), pNS.get());

  auto couplingNSToHeat = CouplingManager::build(COUPLING_TYPE::MEDCOUPLING);
  couplingNSToHeat->setup(pNS.get(), pHeat.get());

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
