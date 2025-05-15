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
  using namespace cocoa;

  std::string configFile = "fd2d_heatpipe.dat";
  if (argc > 1)
  {
    configFile = argv[1];
  }

  auto p = Problem::build(PROBLEM_TYPE::FD2D, EQN_TYPE::HEAT);
  p->setup({{"config_file", configFile}});

  p->print();

  auto const pDer = dynamic_cast<ProblemFD2D *>(p.get());
  auto const dofTop = sideDOF(pDer->mesh_.n_, FD_BC_SIDE::TOP);
  // TODO: do not assume constant profile on inlet
  double const uBottomMean = pDer->bcs_[0].bottom().values[0];

  // // set velocity profile
  // for (uint j = 0U; j < pDer->n_[1]; j++)
  //   for (uint i = 0U; i < pDer->n_[0]; i++)
  //   {
  //     uint const id = i + pDer->n_[0] * j;
  //     double const x = pDer->h_[0] * i;
  //     pDer->c_[1][id] = 1.5 * (1.0 - x * x);
  //   }

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();

    // cyclic bc from top to bottom
    // - get top values
    VectorFD uTop(pDer->mesh_.n_[0]);
    for (uint k = 0U; k < uTop.size(); k++)
    {
      uTop.set(k, pDer->u_[dofTop[k]]);
    }

    // - compute mean top value
    double const uTopMean =
        std::accumulate(uTop.data_.begin(), uTop.data_.end(), 0.0) / uTop.size();

    // - scale back uTop so that its mean is equal to uBottomMean
    VectorFD uBottom(uTop.size());
    for (uint k = 0u; k < uBottom.size(); k++)
    {
      uBottom.set(k, uTop[k] - uTopMean + uBottomMean);
    }

    // - assign new bc values
    pDer->bcs_[0].bottom().values = uBottom;

    fmt::println("deltaTop:    {:.6e}", uTop[uTop.size() - 1] - uTop[0]);
    fmt::println("deltaBottom: {:.6e}", uBottom[uBottom.size() - 1] - uBottom[0]);
    fmt::println("uTop:    {}", uTop);
    fmt::println("uBottom: {}", uBottom);
  }

  return 0;
}
