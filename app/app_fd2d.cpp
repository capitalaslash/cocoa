// local
#include "coupling/coupling_manager.hpp"
#include "problem/problem.hpp"

int main(int argc, char * argv[])
{
  using namespace cocoa;

  std::string configFile = "fd2d_heat.dat";
  if (argc > 1)
  {
    configFile = argv[1];
  }

  auto p = Problem::build(PROBLEM_TYPE::FD2D, EQN_TYPE::HEAT);
  p->setup({{"config_file", configFile}});

  p->print();

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();
  }

  return 0;
}
