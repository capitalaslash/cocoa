// stl
#include <algorithm>

// local
#include "enums.hpp"
#include "la.hpp"
#include "problem/fdutils.hpp"
#include "problem/problem_fd2d.hpp"

auto constexpr N = 80u;
auto constexpr DT = 1e+1;
auto constexpr FINAL_TIME = 40.0;
auto constexpr PRINT_STEP = 1u;

using namespace cocoa;

auto setupFluid(ProblemFD2D * p) -> void
{
  p->name_ = "chtsf";
  p->debug_ = false;
  p->computeCFL_ = true;

  // mesh
  p->mesh_ = MeshFD2D({0.0, 0.0}, {1.0, 0.75}, {N, 3 * N / 4});

  // coupling
  p->couplingType_ = COUPLING_TYPE::MEDCOUPLING;

  // fields
  p->nVars_ = 1u;
  p->varNames_ = {"uf"};
  p->u_.resize(p->mesh_.nPts(), 1.0);
  p->u_ *= 1.0;
  p->uOld_ = p->u_;

  // io setup
  p->outputPrefix_ = "output_chts";
  p->initMeshCoupling();
  p->initFieldCoupling();
  std::filesystem::create_directories(p->outputPrefix_);

  // parameters
  VectorFD cx{p->mesh_.nPts(), 0.0};
  VectorFD cy{p->mesh_.nPts(), 0.0};
  for (uint j = 0u; j < p->mesh_.n_[1]; j++)
    for (uint i = 0u; i < p->mesh_.n_[0]; i++)
    {
      auto const id = i + j * p->mesh_.n_[0];
      auto const x = p->mesh_.pt({i, j});

      cx.set(id, 2.0 * (1.5 * x[1] - 2.0 * x[1] * x[1]));
    }
  p->params_.set("alpha", 0.1);
  p->fields_.emplace("cx", cx);
  p->fields_.emplace("cy", cy);

  p->printFields();

  // bcs
  p->bcs_.resize(p->nVars_);
  p->bcs_[0].left() =
      FDBC{FD_BC_SIDE::LEFT, FD_BC_TYPE::DIRICHLET, 1.0, p->mesh_.n_[1]};
  // p->bcs_[0].left() = FDBC{FD_BC_SIDE::LEFT, FD_BC_TYPE::NEUMANN, 0.0,
  // p->mesh_.n_[1]};
  p->bcs_[0].right() =
      FDBC{FD_BC_SIDE::RIGHT, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[1]};
  p->bcs_[0].bottom() =
      FDBC{FD_BC_SIDE::BOTTOM, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[0]};
  p->bcs_[0].top() = FDBC{FD_BC_SIDE::TOP, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[0]};

  // time
  p->time = 0.0;
  p->finalTime_ = FINAL_TIME;
  p->dt_ = DT;
  p->printStep_ = PRINT_STEP;

  // la
  p->solverType_ = FD_SOLVER_TYPE::BICGSTAB;
  p->maxIters_ = 10000u;
  p->tol_ = 1.0e-6;
  p->m_.init(p->mesh_.nPts(), 5u);
  p->rhs_.resize(p->mesh_.nPts(), 0.0);

  auto const assembly = [](ProblemFD2D * p)
  {
    fmt::println("fluid assembly");

    auto const nx = p->mesh_.n_[0];
    auto const ny = p->mesh_.n_[1];
    auto const iHx = 1.0 / p->mesh_.h_[0];
    auto const iHy = 1.0 / p->mesh_.h_[1];
    auto const iHx2 = iHx * iHx;
    auto const iHy2 = iHy * iHy;
    auto const iDt = p->dt_;

    auto const alpha = p->params_.get<FD_PARAM_TYPE::SCALAR>("alpha");
    auto const cx = p->fields_["cx"];
    auto const cy = p->fields_["cy"];

    for (uint j = 0u; j < ny; j++)
      for (uint i = 0u; i < nx; i++)
      {
        auto const id = i + j * nx;

        // diagonal
        // p->m_.add(id, id, iDt + 2.0 * alpha * (iHx2 + iHy2));
        p->m_.add(id, id, iDt + 2.0 * alpha * (iHx2 + iHy2) + cx[id] * iHx);

        // left
        // auto const valueLeft = -alpha * iHx2 - 0.5 * cx[id] * iHx;
        auto const valueLeft = -alpha * iHx2 - cx[id] * iHx;
        if (i > 0u)
          p->m_.add(id, id - 1, valueLeft);
        else
        {
          p->m_.add(id, id + 1, valueLeft);
          p->bcs_[0].left().ghostValues.set(j, valueLeft);
        }

        // right
        // auto const valueRight = -alpha * iHx2 + 0.5 * cx[id] * iHx;
        auto const valueRight = -alpha * iHx2;
        if (i < nx - 1)
          p->m_.add(id, id + 1, valueRight);
        else
        {
          p->m_.add(id, id - 1, valueRight);
          p->bcs_[0].right().ghostValues.set(j, valueRight);
        }

        // bottom
        auto const valueBottom = -alpha * iHy2 - 0.5 * cy[id] * iHy;
        if (j > 0u)
          p->m_.add(id, id - nx, valueBottom);
        else
        {
          p->m_.add(id, id + nx, valueBottom);
          p->bcs_[0].bottom().ghostValues.set(i, valueBottom);
        }

        // top
        auto const valueTop = -alpha * iHy2 + 0.5 * cy[id] * iHy;
        if (j < ny - 1)
          p->m_.add(id, id + nx, valueTop);
        else
        {
          p->m_.add(id, id - nx, valueTop);
          p->bcs_[0].top().ghostValues.set(i, valueTop);
        }

        // rhs
        p->rhs_.set(id, p->uOld_[id] * iDt);
      }

    p->m_.close();
  };

  p->assemblies_.emplace(EQN_TYPE::CUSTOM, assembly);
  p->eqnType_ = EQN_TYPE::CUSTOM;

  // auto const fixBC = [](ProblemFD2D * p)
  // {
  //   fmt::println("fixing bcs");
  //   for (uint j = 0u; j < p->mesh_.n_[1]; j++)
  //   {
  //     auto const id = 0u + j * p->mesh_.n_[0];
  //     auto const x = p->mesh_.pt({0u, j});

  //     if (x[1] < 0.75)
  //     {
  //       p->m_.clearRow(id);
  //       p->m_.add(id, id, 1.0);
  //       p->rhs_.set(id, 1.0);
  //     }
  //   }
  //   p->m_.close();
  // };
  // p->preSolveFun_.reset(new ProblemFD2D::Assembly_T{fixBC});
}

auto setupSolid(ProblemFD2D * p) -> void
{
  p->name_ = "chtss";
  p->debug_ = false;
  p->computeCFL_ = false;

  // mesh
  p->mesh_ = MeshFD2D({0.0, 0.75}, {1.0, 1.0}, {N, N / 4});

  // coupling
  p->couplingType_ = COUPLING_TYPE::MEDCOUPLING;

  // fields
  p->nVars_ = 1u;
  p->varNames_ = {"us"};
  p->u_.resize(p->mesh_.nPts(), 1.0);
  p->u_ *= 1.0;
  p->uOld_ = p->u_;

  // io setup
  p->outputPrefix_ = "output_chts";
  p->initMeshCoupling();
  p->initFieldCoupling();
  std::filesystem::create_directories(p->outputPrefix_);

  // parameters
  p->params_.set("alpha", 1.0);
  VectorFD q{p->mesh_.nPts(), 6.0};
  p->fields_.emplace("q", q);

  p->printFields();

  // bcs
  p->bcs_.resize(p->nVars_);
  p->bcs_[0].left() = FDBC{FD_BC_SIDE::LEFT, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[1]};
  p->bcs_[0].right() =
      FDBC{FD_BC_SIDE::RIGHT, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[1]};
  p->bcs_[0].bottom() =
      FDBC{FD_BC_SIDE::BOTTOM, FD_BC_TYPE::DIRICHLET, 1.0, p->mesh_.n_[0]};
  p->bcs_[0].top() = FDBC{FD_BC_SIDE::TOP, FD_BC_TYPE::NEUMANN, 0.0, p->mesh_.n_[0]};

  // time
  p->time = 0.0;
  p->finalTime_ = FINAL_TIME;
  p->dt_ = DT;
  p->printStep_ = PRINT_STEP;

  // la
  p->solverType_ = FD_SOLVER_TYPE::BICGSTAB;
  p->maxIters_ = 10000u;
  p->tol_ = 1.0e-6;
  p->m_.init(p->mesh_.nPts(), 5u);
  p->rhs_.resize(p->mesh_.nPts(), 0.0);

  auto const assembly = [](ProblemFD2D * p)
  {
    fmt::println("solid assembly");

    auto const nx = p->mesh_.n_[0];
    auto const ny = p->mesh_.n_[1];
    auto const iHx = 1.0 / p->mesh_.h_[0];
    auto const iHy = 1.0 / p->mesh_.h_[1];
    auto const iHx2 = iHx * iHx;
    auto const iHy2 = iHy * iHy;
    auto const iDt = p->dt_;

    auto const alpha = p->params_.get<FD_PARAM_TYPE::SCALAR>("alpha");
    auto const q = p->fields_["q"];

    for (uint j = 0u; j < ny; j++)
      for (uint i = 0u; i < nx; i++)
      {
        auto const id = i + j * nx;

        // diagonal
        p->m_.add(id, id, iDt + 2.0 * alpha * (iHx2 + iHy2));

        // left
        auto const valueLeft = -alpha * iHx2;
        if (i > 0)
          p->m_.add(id, id - 1, valueLeft);
        else
        {
          p->m_.add(id, id + 1, valueLeft);
          p->bcs_[0].left().ghostValues.set(j, valueLeft);
        }

        // right
        auto const valueRight = -alpha * iHx2;
        if (i < nx - 1)
          p->m_.add(id, id + 1, valueRight);
        else
        {
          p->m_.add(id, id - 1, valueRight);
          p->bcs_[0].right().ghostValues.set(j, valueRight);
        }

        // bottom
        auto const valueBottom = -alpha * iHy2;
        if (j > 0u)
          p->m_.add(id, id - nx, valueBottom);
        else
        {
          p->m_.add(id, id + nx, valueBottom);
          p->bcs_[0].bottom().ghostValues.set(i, valueBottom);
        }

        // top
        auto const valueTop = -alpha * iHy2;
        if (j < ny - 1)
          p->m_.add(id, id + nx, valueTop);
        else
        {
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
}

auto getSolOnSide(ProblemFD2D const * p, FD_BC_SIDE const side) -> VectorFD
{
  auto const dofs = sideDOF(p->mesh_.n_, side);
  auto solSide = VectorFD{dofs.size()};
  for (auto i = 0u; i < dofs.size(); i++)
    solSide.set(i, p->u_[dofs[i]]);
  return solSide;
}

auto getGradOnSide(ProblemFD2D const * p, FD_BC_SIDE const side) -> VectorFD
{
  auto const dofs = sideDOF(p->mesh_.n_, side);
  auto const offset = sideOffset(p->mesh_.n_, side);
  auto const h = (side == FD_BC_SIDE::LEFT || side == FD_BC_SIDE::RIGHT)
                     ? p->mesh_.h_[0]
                     : p->mesh_.h_[1];
  auto gradSide = VectorFD{dofs.size()};
  for (auto i = 0u; i < dofs.size(); i++)
  {
    // // linear extrapolation
    auto const valueL = (p->u_[dofs[i]] - p->u_[dofs[i] + offset]) / h;
    // quadratic extrapolation
    auto const valueQ = (1.5 * p->u_[dofs[i]] - 2. * p->u_[dofs[i] + offset] +
                         0.5 * p->u_[dofs[1] + 2 * offset]) /
                        h;

    auto const value = (valueQ / valueL > 2) ? valueL : valueQ;

    gradSide.set(i, value);
  }
  return gradSide;
}

int main()
{
  auto pFluid = std::unique_ptr<ProblemFD2D>{new ProblemFD2D};
  setupFluid(pFluid.get());
  pFluid->print();

  auto pSolid = std::unique_ptr<ProblemFD2D>{new ProblemFD2D};
  setupSolid(pSolid.get());
  pSolid->print();

  VectorFD uI{pSolid->mesh_.n_[0], 1.0};
  VectorFD fluxI{pSolid->mesh_.n_[0], -15.0};
  // VectorFD uFixed{pSolid->mesh_.n_[0], 1.0};
  // for (uint k = 0u; k < pSolid->mesh_.n_[0]; k++)
  //   uFixed.set(k, 2. + 2. * k / N);
  // VectorFD gFixed{pSolid->mesh_.n_[0], -15.};
  while (pFluid->run() && pSolid->run())
  {
    pSolid->advance();
    pFluid->advance();

    for (uint i = 0u; i < 3u; i++)
    {
      // fluid
      pFluid->bcs_[0].top().values = fluxI;
      // pFluid->bcs_[0].top().values = gFixed;

      auto numItersFluid = 100u;
      for (auto k = 0u; k < 100u; k++)
      {
        numItersFluid = pFluid->solve();
        if (std::fabs(pFluid->u_.norm2sq() - pFluid->uOld_.norm2sq()) < 1e-4)
        {
          fmt::print("fluid internal iterations: {}\n", k);
          break;
        }
      }

      uI = getSolOnSide(pFluid.get(), FD_BC_SIDE::TOP);
      auto const uMin = std::min_element(uI.data_.begin(), uI.data_.end());
      auto const uMax = std::max_element(uI.data_.begin(), uI.data_.end());
      fmt::print("uI min: {:.6e}\n", *uMin);
      fmt::print("uI max: {:.6e}\n", *uMax);

      pSolid->bcs_[0].bottom().values = uI;
      // pSolid->bcs_[0].bottom().values = uFixed;

      auto numItersSolid = 100u;
      for (auto k = 0u; k < 100u; k++)
      {
        numItersSolid = pSolid->solve();
        if (std::fabs(pSolid->u_.norm2sq() - pSolid->uOld_.norm2sq()) < 1e-4)
        {
          fmt::print("internal solid iterations: {}\n", k);
          break;
        }
      }

      fluxI = getGradOnSide(pSolid.get(), FD_BC_SIDE::BOTTOM);
      fluxI *= pFluid->params_.get<FD_PARAM_TYPE::SCALAR>("alpha") /
               pSolid->params_.get<FD_PARAM_TYPE::SCALAR>("alpha");
      auto const fluxMin = std::min_element(fluxI.data_.begin(), fluxI.data_.end());
      auto const fluxMax = std::max_element(fluxI.data_.begin(), fluxI.data_.end());
      fmt::print("fluxI min: {:.6e}\n", *fluxMin);
      fmt::print("fluxI max: {:.6e}\n", *fluxMax);

      auto totalFlux = 0.0;
      totalFlux += fluxI[0] * 0.5 * pSolid->mesh_.h_[0];
      for (uint k = 1u; k < pSolid->mesh_.n_[0] - 1; k++)
        totalFlux += fluxI[k] * pSolid->mesh_.h_[0];
      totalFlux += fluxI[pSolid->mesh_.n_[0] - 1] * 0.5 * pSolid->mesh_.h_[0];
      fmt::print("totalFlux: {:.6e}\n", totalFlux);

      // if (numItersSolid + numItersFluid == 2u)
      // {
      //   fmt::print("breaking after {} internal iterations\n", i);
      //   break;
      // }
    }

    pSolid->print();
    pFluid->print();
  }

  return 0;
}
