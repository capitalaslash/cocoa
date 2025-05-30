// local
#include "problem/problem.hpp"

int main()
{
  using namespace cocoa;

  auto p = Problem::build(PROBLEM_TYPE::OFORG, EQN_TYPE::NONE);
  p->setup({
      {"case_dir", std::filesystem::path{"oforg1"}},
      {"config_file", std::filesystem::path{"oforg1.dat"}},
  });

  p->print();

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();
  }

  return 0;
}
