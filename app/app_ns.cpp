// local
#include "problem/problem.hpp"

int main()
{
  using namespace cocoa;

  auto p = Problem::build(PROBLEM_TYPE::PROXPDE, EQN_TYPE::NS);
  p->setup({{"config_file", std::filesystem::path{"proxpde_ns.yaml"}}});

  p->print();

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();
  }

  return 0;
}
