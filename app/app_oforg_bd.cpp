// fmt
#include <fmt/ranges.h>
#include <fmt/std.h>

// local
#include "coupling/coupling_manager.hpp"
#include "coupling/coupling_med.hpp"
#include "problem/problem.hpp"

int main()
{
  using namespace cocoa;

  auto p1 = Problem::build(PROBLEM_TYPE::OFORG, EQN_TYPE::NONE);
  p1->setup({
      {"case_dir", std::filesystem::path{"oforg_channel1"}},
      {"config_file", std::filesystem::path{"oforg_channel1.dat"}},
  });

  auto p2 = Problem::build(PROBLEM_TYPE::OFORG, EQN_TYPE::NONE);
  p2->setup({
      {"case_dir", std::filesystem::path{"oforg_channel2"}},
      {"config_file", std::filesystem::path{"oforg_channel2.dat"}},
  });

  p1->print();
  p2->print();

  auto coupling12 =
      CouplingManager::build(COUPLING_TYPE::MEDCOUPLING, COUPLING_SCOPE::BOUNDARY);
  coupling12->setup(
      {p1.get(), "top", {"U"}},
      {p2.get(), "bottom", {"U"}},
      INTERPOLATION_METHOD::P0P0);

  auto coupling21 =
      CouplingManager::build(COUPLING_TYPE::MEDCOUPLING, COUPLING_SCOPE::BOUNDARY);
  coupling21->setup(
      {p2.get(), "bottom", {"p"}},
      {p1.get(), "top", {"p"}},
      INTERPOLATION_METHOD::P0P0);

  while (p1->run() || p2->run())
  {
    p1->advance();
    p2->advance();

    if (p1->time > 1e2)
      coupling21->project("p");
    p1->solve();
    if (p2->time > 1e2)
      coupling12->project("U");
    p2->solve();

    p1->print();
    p2->print();
  }

  return 0;
}
