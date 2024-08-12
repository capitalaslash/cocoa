#include <vector>

#include "coupling_manager.hpp"
#include "problem_factory.hpp"

int main()
{

  // ProblemOForg p1;
  // manager.problems_.push_back(&p1);

  // ProblemFEMUS p2;
  // manager.problems_.push_back(&p2);

  // ProblemProXPDE p1;
  // ProblemProXPDE p2;

  std::unique_ptr<Problem> p1{buildProblem(PROBLEM_TYPE::FD1D)};
  std::unique_ptr<Problem> p2{buildProblem(PROBLEM_TYPE::FD1D)};

  CouplingManager coupling12;
  coupling12.setup(p1.get(), p2.get());
  CouplingManager coupling21;
  coupling21.setup(p2.get(), p1.get());

  // p1->setup({});
  // p2->setup({});
  // p1->setup({{"config_file", "proxpde1.yaml"}});
  // p2->setup({{"config_file", "proxpde2.yaml"}});
  p1->setup({{"config_file", "fd1d1.dat"}});
  p2->setup({{"config_file", "fd1d2.dat"}});

  p1->print();
  p2->print();

  while (p1->run() || p2->run())
  {
    // coupling21.project("Tsys", "T");
    p1->advance();
    p1->solve();
    p1->print();

    coupling12.project("T", "Tsys");
    p2->advance();
    p2->solve();
    p2->print();
  }

  return 0;
}
