// local
#include "problem/problem.hpp"

int main()
{
  auto p = Problem::build(PROBLEM_TYPE::FD1D, EQN_TYPE::HEAT);
  p->setup({{"config_file", "fd1d_heat.dat"}});

  p->print();

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();
  }

  auto const * sol = p->getField("T");
  auto const computed = sol->at(sol->size() - 1);
  auto const expected = 2.3956041663411976e-06;
  if (std::fabs(computed - expected) > 1.e-12)
  {
    fmt::print(
        stderr,
        "computed value ({:.16e}) is different from the expected {:.16e}\n",
        computed,
        expected);
    return 1;
  }
  return 0;
}
