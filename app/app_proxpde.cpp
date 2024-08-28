// local
#include "coupling_factory.hpp"
#include "problem_factory.hpp"

int main()
{
  std::unique_ptr<Problem> p1{buildProblem(PROBLEM_TYPE::PROXPDE)};
  p1->setup({{"config_file", "proxpde1.yaml"}});
  dynamic_cast<ProblemProXPDE *>(p1.get())->q_
      << [](proxpde::Vec3 const & p) { return std::sin(M_PI * p[1]); };

  std::unique_ptr<Problem> p2{buildProblem("proxpde")};
  p2->setup({{"config_file", "proxpde2.yaml"}});

  // check that both problems have been set to use the same coupling type
  assert(p1->couplingType_ == p2->couplingType_);

  std::unique_ptr<CouplingManager> coupling12{buildCoupling(p1->couplingType_)};
  coupling12->setup(p1.get(), p2.get());

  p1->print();
  p2->print();

  while (p1->run() || p2->run())
  {
    p1->advance();
    p1->solve();
    p1->print();

    coupling12->project("T", "Tcfd");
    p2->advance();
    p2->solve();
    p2->print();
  }

  return 0;
}
