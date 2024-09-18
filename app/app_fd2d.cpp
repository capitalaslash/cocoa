// local
#include "problem/problem.hpp"

int main()
{
  auto p = Problem::build(PROBLEM_TYPE::FD2D, EQN_TYPE::HEAT);
  p->setup({{"config_file", "fd2d_heat.dat"}});

  p->print();

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();
  }

  return 0;
}
