// local
#include "enums.hpp"
#include "problem/problem_fd2d.hpp"

constexpr auto U_RIGHT = 0.0;

using namespace cocoa;

auto setup(ProblemFD2D * p) -> void
{
  p->name_ = "fd2d_custom";
  p->debug_ = false;

  // mesh
  p->mesh_ = MeshFD2D({0.0, 0.0}, {1.0, 1.0}, {20u, 20u});

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
  p->bcs_[0].left() =
      FDBC{FD_BC_SIDE::LEFT, FD_BC_TYPE::DIRICHLET, 0.0, p->mesh_.n_[1]};
  p->bcs_[0].right() =
      FDBC{FD_BC_SIDE::RIGHT, FD_BC_TYPE::NEUMANN, 0.5 - U_RIGHT, p->mesh_.n_[1]};
  p->bcs_[0].top() = FDBC{FD_BC_SIDE::TOP, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[0]};
  p->bcs_[0].bottom() =
      FDBC{FD_BC_SIDE::BOTTOM, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[0]};

  // time
  p->time = 0.0;
  p->finalTime_ = 20.0;
  p->dt_ = 1.0;

  // la
  p->solverType_ = FD_SOLVER_TYPE::CG;
  p->maxIters_ = 10000u;
  p->tol_ = 1.0e-6;
  p->m_.init(p->mesh_.nPts(), 5u);
  p->rhs_.resize(p->mesh_.nPts(), 0.0);

  auto const assembly = [](ProblemFD2D * p)
  {
    fmt::println("custom assembly");

    auto const iDt = 1.0 / p->dt_;
    auto const nx = p->mesh_.n_[0];
    auto const ny = p->mesh_.n_[1];
    auto const iH2x = 1.0 / (p->mesh_.h_[0] * p->mesh_.h_[0]);
    auto const iH2y = 1.0 / (p->mesh_.h_[1] * p->mesh_.h_[1]);
    auto const alpha = p->params_.get<FD_PARAM_TYPE::SCALAR>("alpha");
    auto const q = p->fields_["q"];

    for (uint j = 0u; j < ny; j++)
      for (uint i = 0u; i < nx; i++)
      {
        auto const id = i + j * nx;

        // middle
        double const value = iDt + 2. * alpha * (iH2x + iH2y);
        p->m_.add(id, id, value);

        // bottom
        double const valueBottom = -alpha * iH2y;
        if (j > 0u)
          p->m_.add(id, id - nx, valueBottom);
        else
        {
          p->m_.add(id, id + nx, valueBottom);
          p->bcs_[0].bottom().ghostValues.set(i, valueBottom);
        }

        // right
        double const valueRight = -alpha * iH2x;
        if (i < nx - 1)
          p->m_.add(id, id + 1, valueRight);
        else
        {
          p->m_.add(id, id - 1, valueRight);
          p->bcs_[0].right().ghostValues.set(j, valueRight);
        }

        // top
        double const valueTop = -alpha * iH2y;
        if (j < ny - 1)
          p->m_.add(id, id + nx, valueTop);
        else
        {
          p->m_.add(id, id - nx, valueTop);
          p->bcs_[0].top().ghostValues.set(i, valueTop);
        }

        // left
        double const valueLeft = -alpha * iH2x;
        if (i > 0u)
          p->m_.add(id, id - 1, valueLeft);
        else
        {
          p->m_.add(id, id + 1, valueLeft);
          p->bcs_[0].left().ghostValues.set(j, valueLeft);
        }

        // rhs
        p->rhs_.set(id, p->uOld_[id] * iDt + q[id]);
      }

    p->m_.close();
  };

  p->assemblies_.emplace(EQN_TYPE::CUSTOM, assembly);
  p->eqnType_ = EQN_TYPE::CUSTOM;

  p->outputPrefix_ = "output_fd2d_custom";
  p->cleanOutput_ = true;
  p->initOutput();
}

int main()
{
  auto p = std::unique_ptr<ProblemFD2D>{new ProblemFD2D};
  setup(p.get());

  p->print();

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();
  }

  auto pDer = dynamic_cast<ProblemFD2D *>(p.get());
  auto const & sol = pDer->u_;
  auto const computed = sol[sol.size() - 1];
  auto const expected = -3.3533540879997921e-09 + U_RIGHT;
  if (std::fabs(computed - expected) > 1.0e-12)
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
