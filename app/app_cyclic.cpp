// std
#include <numeric>

// fmt
#include <fmt/ranges.h>

// local
#include "problem/problem.hpp"

// tmp
#include "problem/problem_fd2d.hpp"

int main(int argc, char * argv[])
{
  std::string configFile = "fd2d_heatpipe.dat";
  if (argc > 1)
  {
    configFile = argv[1];
  }

  auto p = Problem::build(PROBLEM_TYPE::FD2D, EQN_TYPE::HEAT);
  p->setup({{"config_file", configFile}});

  p->print();

  auto const pDer = dynamic_cast<ProblemFD2D *>(p.get());
  auto const dofTop = sideDOF(pDer->n_, 2 /*TOP*/);
  // TODO: do not assume constant profile on inlet
  double const uBottomMean = pDer->bcs_[0].values[0];

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();

    // cyclic bc from top to bottom
    // - get top values
    std::vector<double> uTop(pDer->n_[0]);
    for (uint k = 0U; k < uTop.size(); k++)
    {
      uTop[k] = pDer->u_[dofTop[k]];
    }

    // - compute mean top value
    double const uTopMean =
        std::accumulate(uTop.begin(), uTop.end(), 0.0) / uTop.size();

    // - scale back uTop so that its mean is equal to uBottomMean
    std::vector<double> uBottom(uTop);
    for (auto & value: uBottom)
    {
      value -= uTopMean + uBottomMean;
    }

    // - assign new bc values
    pDer->bcs_[0].values = uBottom;

    fmt::print("deltaTop:    {:.6e}\n", uTop[uTop.size() - 1] - uTop[0]);
    fmt::print("deltaBottom: {:.6e}\n", uBottom[uBottom.size() - 1] - uBottom[0]);
    fmt::print("uTop:    {::+.6e}\n", uTop);
    fmt::print("uBottom: {::+.6e}\n", uBottom);
  }

  return 0;
}
