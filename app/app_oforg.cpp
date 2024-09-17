// local
#include "problem/problem.hpp"

int main()
{
  auto p = Problem::build(PROBLEM_TYPE::OFORG, EQN_TYPE::HEAT);
  p->setup({{"case_dir", "oforg1"}});

  p->print();

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();
  }

  return 0;
}
