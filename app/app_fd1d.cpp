// local
#include "problem/problem.hpp"

int main()
{
  auto p = Problem::build(PROBLEM_TYPE::FD1D, EQN_TYPE::HEAT);
  p->setup({{"config_file", "fd1d1.dat"}});

  p->print();

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();
  }

  return 0;
}
