// local
#include "problem/problem.hpp"
#include "problem/problem_fd1d.hpp"

int main()
{
  using namespace cocoa;

  auto p = Problem::build(PROBLEM_TYPE::FD1D, EQN_TYPE::HEAT);
  p->setup({{"config_file", std::filesystem::path{"fd1d_heat.dat"}}});

  p->print();

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();
  }

  auto pDer = dynamic_cast<ProblemFD1D *>(p.get());
  auto const & sol = pDer->u_;
  auto const computed = sol[sol.size() - 1];
  auto const expected = -2.5330183695843534e-06;
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
