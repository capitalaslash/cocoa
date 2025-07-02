// local
#include "enums.hpp"
#include "problem/fdutils.hpp"
#include "problem/problem_fd2d.hpp"

// // arithmetic
// static constexpr auto mean(double const a, double const b) -> double { return 0.5 *
// (a + b); }

// harmonic
static constexpr auto mean(double const a, double const b) -> double
{
  return 2 * a * b / (a + b);
}

using namespace cocoa;

auto setup(ProblemFD2D * p) -> void
{
  p->name_ = "chtm";
  p->debug_ = false;
  p->computeCFL_ = true;

  // mesh
  p->mesh_ = MeshFD2D({0.0, 0.0}, {1.0, 1.0}, {80u, 80u});

  // fields
  p->nVars_ = 1u;
  p->varNames_ = {"u"};
  p->u_.resize(p->mesh_.nPts(), 1.0);
  p->u_ *= 1.0;
  p->uOld_ = p->u_;

  // io setup
  p->outputPrefix_ = "output_chtm";
  std::filesystem::create_directories(p->outputPrefix_);

  // parameters
  VectorFD alpha{p->mesh_.nPts(), 0.0};
  VectorFD q{p->mesh_.nPts(), 0.0};
  VectorFD cx{p->mesh_.nPts(), 0.0};
  VectorFD cy{p->mesh_.nPts(), 0.0};
  auto const alphaFluid = 0.1;
  auto const alphaSolid = 1.0;
  auto const interfaceY = 0.75;
  for (uint j = 0u; j < p->mesh_.n_[1]; j++)
    for (uint i = 0u; i < p->mesh_.n_[0]; i++)
    {
      auto const id = i + j * p->mesh_.n_[0];
      auto const x = p->mesh_.pt({i, j});

      // fluid
      if (x[1] < interfaceY - 1e-6)
      {
        q.set(id, 0.0);
        alpha.set(id, alphaFluid);
        cx.set(id, 2.0 * (1.5 * x[1] - 2.0 * x[1] * x[1]));
      }
      // solid
      else if (x[1] > interfaceY + 1e-6)
      {
        q.set(id, 6.0);
        alpha.set(id, alphaSolid);
      }
      // interface
      else
      {
        q.set(id, 3.0);
        alpha.set(id, mean(alphaFluid, alphaSolid));
      }
    }
  p->fields_.emplace("alpha", alpha);
  p->fields_.emplace("q", q);
  p->fields_.emplace("cx", cx);
  p->fields_.emplace("cy", cy);

  p->printFields();

  // bcs
  p->bcs_.resize(p->nVars_);
  p->bcs_[0].left() = FDBC{FD_BC_SIDE::LEFT, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[1]};
  p->bcs_[0].right() =
      FDBC{FD_BC_SIDE::RIGHT, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[1]};
  p->bcs_[0].bottom() =
      FDBC{FD_BC_SIDE::BOTTOM, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[0]};
  p->bcs_[0].top() = FDBC{FD_BC_SIDE::TOP, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[0]};

  // time
  p->time = 0.0;
  p->finalTime_ = 20.0;
  p->dt_ = 1.0;
  p->printStep_ = 1u;

  // la
  p->solverType_ = FD_SOLVER_TYPE::BICGSTAB;
  p->maxIters_ = 10000u;
  p->tol_ = 1.0e-6;
  p->m_.init(p->mesh_.nPts(), 5u);
  p->rhs_.resize(p->mesh_.nPts(), 0.0);

  auto const assembly = [](ProblemFD2D * p)
  {
    fmt::println("cht assembly");

    auto const nx = p->mesh_.n_[0];
    auto const ny = p->mesh_.n_[1];
    auto const iHx = 1.0 / p->mesh_.h_[0];
    auto const iHy = 1.0 / p->mesh_.h_[1];
    auto const iHx2 = iHx * iHx;
    auto const iHy2 = iHy * iHy;
    auto const iDt = p->dt_;

    auto const alpha = p->fields_["alpha"];
    auto const q = p->fields_["q"];
    auto const cx = p->fields_["cx"];
    auto const cy = p->fields_["cy"];

    for (uint j = 0u; j < ny; j++)
      for (uint i = 0u; i < nx; i++)
      {
        auto const id = i + j * nx;

        // diagonal
        if (i > 0u && i < nx - 1 && j > 0u && j < ny - 1)
        {
          auto const value =
              iDt +
              (mean(alpha[id - 1], alpha[id]) + mean(alpha[id + 1], alpha[id])) * iHx2 +
              (mean(alpha[id - nx], alpha[id]) + mean(alpha[id + nx], alpha[id])) *
                  iHy2;
          p->m_.add(id, id, value);
        }
        else
          p->m_.add(id, id, iDt + 2.0 * alpha[id] * (iHx2 + iHy2));

        // left
        if (i > 0u && i < nx - 1 && j > 0u && j < ny - 1)
        { // internal
          auto const valueLeft =
              -mean(alpha[id - 1], alpha[id]) * iHx2 - 0.5 * cx[id] * iHx;
          p->m_.add(id, id - 1, valueLeft);
        }
        else if (i > 0u)
        { // other sides
          auto const valueLeft = -alpha[id] * iHx2 - 0.5 * cx[id] * iHx;
          p->m_.add(id, id - 1, valueLeft);
        }
        else
        { // own side
          auto const valueLeft = -alpha[id] * iHx2 - 0.5 * cx[id] * iHx;
          p->m_.add(id, id + 1, valueLeft);
          p->bcs_[0].left().ghostValues.set(j, valueLeft);
        }

        // right
        if (i > 0u && i < nx - 1 && j > 0u && j < ny - 1)
        { // internal
          auto const valueRight =
              -mean(alpha[id + 1], alpha[id]) * iHx2 + 0.5 * cx[id] * iHx;
          p->m_.add(id, id + 1, valueRight);
        }
        else if (i < nx - 1)
        { // other sides
          auto const valueRight = -alpha[id] * iHx2 + 0.5 * cx[id] * iHx;
          p->m_.add(id, id + 1, valueRight);
        }
        else
        { // own side
          auto const valueRight = -alpha[id] * iHx2 + 0.5 * cx[id] * iHx;
          p->m_.add(id, id - 1, valueRight);
          p->bcs_[0].right().ghostValues.set(j, valueRight);
        }

        // bottom
        if (i > 0u && i < nx - 1 && j > 0u && j < ny - 1)
        { // internal
          auto const valueBottom =
              -mean(alpha[id - nx], alpha[id]) * iHy2 - 0.5 * cy[id] * iHy;
          p->m_.add(id, id - nx, valueBottom);
        }
        else if (j > 0u)
        { // other sides
          auto const valueBottom = -alpha[id] * iHy2 - 0.5 * cy[id] * iHy;
          p->m_.add(id, id - nx, valueBottom);
        }
        else
        { // own side
          auto const valueBottom = -alpha[id] * iHy2 - 0.5 * cy[id] * iHy;
          p->m_.add(id, id + nx, valueBottom);
          p->bcs_[0].bottom().ghostValues.set(i, valueBottom);
        }

        // top
        if (i > 0u && i < nx - 1 && j > 0u && j < ny - 1)
        { // internal
          auto const valueTop =
              -mean(alpha[id + nx], alpha[id]) * iHy2 + 0.5 * cy[id] * iHy;
          p->m_.add(id, id + nx, valueTop);
        }
        else if (j < ny - 1)
        { // other sides
          auto const valueTop = -alpha[id] * iHy2 + 0.5 * cy[id] * iHy;
          p->m_.add(id, id + nx, valueTop);
        }
        else
        { // own side
          auto const valueTop = -alpha[id] * iHy2 + 0.5 * cy[id] * iHy;
          p->m_.add(id, id - nx, valueTop);
          p->bcs_[0].top().ghostValues.set(i, valueTop);
        }

        // rhs
        p->rhs_.set(id, p->uOld_[id] * iDt + q[id]);
      }

    p->m_.close();
  };

  p->assemblies_.emplace(EQN_TYPE::CUSTOM, assembly);
  p->eqnType_ = EQN_TYPE::CUSTOM;

  auto const fixBC = [](ProblemFD2D * p)
  {
    fmt::println("fixing bcs");
    for (uint j = 0u; j < p->mesh_.n_[1]; j++)
    {
      auto const id = 0u + j * p->mesh_.n_[0];
      auto const x = p->mesh_.pt({0u, j});

      if (x[1] <= 0.75)
      {
        p->m_.clearRow(id);
        p->m_.add(id, id, 1.0);
        p->rhs_.set(id, 1.0);
      }
    }
    p->m_.close();
  };
  p->preSolveFun_.reset(new ProblemFD2D::Assembly_T{fixBC});
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

  return 0;
}
