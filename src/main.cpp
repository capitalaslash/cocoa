#include <vector>

#include <fmt/core.h>
#include <fmt/std.h>

#include "med_manager.hpp"
#include "problem_fd1d.hpp"
#include "problem_femus.hpp"
#include "problem_oforg.hpp"
#include "problem_proxpde.hpp"

int main()
{

  // ProblemOForg p1;
  // manager.problems_.push_back(&p1);

  // ProblemFEMUS p2;
  // manager.problems_.push_back(&p2);

  ProblemProXPDE p3;
  ProblemProXPDE p4;

  ProblemFD1D p5;

  MEDManager coupling34{&p3, &p4};
  MEDManager coupling43{&p4, &p3};

  // p1.setup({});
  // p2.setup({});

  p3.setup({{"config_file", "proxpde1.yaml"}});
  p4.setup({{"config_file", "proxpde2.yaml"}});
  p5.setup({{"config_file", "fd1d.dat"}});

  p3.print();
  p4.print();
  p5.print();

  while (p3.run() || p4.run() || p5.run())
  {
    // coupling43.project("Tsys", "T");
    p3.advance();
    p3.solve();
    p3.print();

    // coupling34.project("T", "Tsys");
    p4.advance();
    p4.solve();
    p4.print();

    p5.advance();
    p5.solve();
    p5.print();
  }

  return 0;
}
