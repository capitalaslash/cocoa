// local
#include "enums.hpp"
#include "problem/problem_fd1d.hpp"

constexpr auto U_RIGHT = 0.0;

using namespace cocoa;

auto setup(ProblemFD1D * p) -> void
{
  p->name_ = "fd1d_custom";
  p->debug_ = false;

  // mesh
  p->mesh_ = MeshFD1D({0.0}, {1.0}, {100});

  // coupling
  p->couplingType_ = COUPLING_TYPE::MEDCOUPLING;

  // fields
  p->nVars_ = 1u;
  p->varNames_ = {"u"};
  p->u_.resize(p->mesh_.nPts(), 1.0);
  p->u_ *= 0.0;
  p->uOld_ = p->u_;

  // eqn params
  p->params_.set("alpha", 1.0);
  VectorFD q(p->mesh_.nPts(), 1.0);
  q *= 1.0;
  p->fields_.emplace("q", q);

  // bcs
  p->bcs_.resize(p->nVars_);
  p->bcs_[0].left() = FDBC{FD_BC_SIDE::LEFT, FD_BC_TYPE::DIRICHLET, 0.0};
  p->bcs_[0].right() = FDBC{FD_BC_SIDE::RIGHT, FD_BC_TYPE::NEUMANN, U_RIGHT - 0.5};

  // time
  p->time = 0.0;
  p->finalTime_ = 20.0;
  p->dt_ = 1.0;

  // la
  p->solverType_ = FD_SOLVER_TYPE::CG;
  p->maxIters_ = 10000u;
  p->tol_ = 1.0e-6;
  p->m_.init(p->mesh_.nPts(), 3u);
  p->rhs_.resize(p->mesh_.nPts(), 0.0);

  auto const assembly = [](ProblemFD1D * p)
  {
    fmt::println("custom assembly");

    auto const iH2 = 1.0 / (p->mesh_.h_[0] * p->mesh_.h_[0]);
    auto const alpha = p->params_.get<FD_PARAM_TYPE::SCALAR>("alpha");
    auto const q = p->fields_["q"];

    for (auto k = 0u; k < p->mesh_.nPts(); k++)
    {
      // diagonal
      p->m_.add(k, k, 1.0 + 2.0 * alpha * p->dt_ * iH2);

      // left
      auto const valueLeft = -alpha * p->dt_ * iH2;
      if (k > 0)
        p->m_.add(k, k - 1, valueLeft);
      else
      {
        p->m_.add(k, k + 1, valueLeft);
        p->bcs_[0].left().ghostValues.set(0, valueLeft);
      }

      // right
      auto const valueRight = -alpha * p->dt_ * iH2;
      if (k < p->mesh_.nPts() - 1)
        p->m_.add(k, k + 1, valueRight);
      else
      {
        p->m_.add(k, k - 1, valueRight);
        p->bcs_[0].right().ghostValues.set(0, valueRight);
      }

      // rhs
      p->rhs_.set(k, p->uOld_[k] + q[k] * p->dt_);
    }

    p->m_.close();
  };

  p->assemblies_.emplace(EQN_TYPE::CUSTOM, assembly);
  p->eqnType_ = EQN_TYPE::CUSTOM;

  p->outputPrefix_ = "output_fd1d_custom";
  p->initMeshCoupling();
  p->initFieldCoupling();
  std::filesystem::create_directories(p->outputPrefix_);
}

int main()
{
  auto p = std::unique_ptr<ProblemFD1D>{new ProblemFD1D};
  // p->setup({{"config_file", std::filesystem::path{"fd1d_heat.dat"}}});
  setup(p.get());

  p->print();

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();
  }

  auto const * sol = p->getField("u");
  auto const computed = sol->at(sol->size() - 1);
  auto const expected = U_RIGHT - 8.1645061956362504e-09;
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
