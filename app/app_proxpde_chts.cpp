// stl
#include <unordered_map>

// libfmt
#include <fmt/ranges.h>

// local
#include "problem/problem_proxpde.hpp"

struct AssemblySolid: public ProblemProXPDE::Assembly
{
  using Elem_T = ProblemProXPDEHeat::Elem_T;
  using FESpaceGradElem_T = proxpde::FESpace<
      ProblemProXPDEHeat::Mesh_T,
      ProblemProXPDEHeat::FE0_T::RefFE_T,
      ProblemProXPDEHeat::FE_T::RecommendedQR,
      Elem_T::dim>;
  static constexpr auto numPtsQRBC =
      proxpde::SideQR_T<typename ProblemProXPDEHeat::FESpace_T::QR_T>::numPts;
  using FESpaceGradBC_T = proxpde::FESpace<
      ProblemProXPDEHeat::Mesh_T,
      ProblemProXPDEHeat::FE0_T::RefFE_T,
      proxpde::SideGaussQR<ProblemProXPDEHeat::Elem_T, numPtsQRBC>,
      Elem_T::dim>;

  auto evaluate(ProblemProXPDE * pParent, proxpde::Builder<> & b) -> void override
  {
    fmt::print("cht solid assembly\n");
    auto p = dynamic_cast<ProblemProXPDEHeat *>(pParent);

    auto const alpha = p->params_["alpha"];
    auto & q = p->fieldsP0_.at("q");

    // assembly lhs
    b.buildLhs(
        std::tuple{
            proxpde::AssemblyMass{1. / p->dt_, p->feSpace_},
            proxpde::AssemblyStiffness{alpha, p->feSpace_},
        },
        p->bcs_);
    b.closeMatrix();

    // assembly rhs
    b.buildRhs(
        std::tuple{
            proxpde::AssemblyProjection{1. / p->dt_, p->uOld_, p->feSpace_},
            proxpde::AssemblyProjection{1., q.data, *q.feSpace, p->feSpace_},
        },
        p->bcs_);
  }

  FESpaceGradElem_T feSpaceGrad_;
  proxpde::FEVar<FESpaceGradElem_T> uGrad_;
  FESpaceGradBC_T feSpaceGradBC_;
  proxpde::FEVar<FESpaceGradBC_T> uGradBC_;
};

auto setupSolid(ProblemProXPDEHeat * p, proxpde::ParameterDict & config) -> void
{
  // name
  p->name_ = "proxpde_chts";

  // parameters
  auto const n = config["n"].as<uint>();
  config["mesh"]["type"] = proxpde::MeshType::STRUCTURED;
  config["mesh"]["origin"] = proxpde::Vec3{0.0, 0.75, 0.0};
  config["mesh"]["length"] = proxpde::Vec3{1.0, 0.25, 0.0};
  config["mesh"]["n"] = std::array{n, n / 4, 0u};
  config["mesh"]["flags"] = proxpde::Bitmask{proxpde::MeshFlags::NORMALS};

  // mesh
  proxpde::readMesh(p->mesh_, config["mesh"]);

  p->initMeshMED(p->name_, p->mesh_);

  // coupling
  p->couplingType_ = COUPLING_TYPE::MEDCOUPLING;

  // time
  p->time = 0.0;
  p->dt_ = config["dt"].as<double>();
  p->finalTime_ = config["final_time"].as<double>();
  p->printStep_ = config["print_step"].as<uint>();

  // finite elements
  p->feSpace_.init(p->mesh_);
  p->feSpaceP0_.init(p->mesh_);
  p->T_.init("T", p->feSpace_);
  p->T_ << 2.0;
  p->u_ = p->T_.data;

  // fields
  p->params_["alpha"] = 1.0;
  p->fieldsP0_.emplace("q", proxpde::FEVar{"q", p->feSpaceP0_});
  p->fieldsP0_.at("q") << 6.0;

  // bcs
  p->bcs_ = std::vector{
      proxpde::BCEss{
          p->feSpace_,
          proxpde::side::BOTTOM,
          [](proxpde::Vec3 const &) { return 2.0; }},
  };

  // io
  p->io_.init(p->feSpace_, "output_" + p->name_ + "/Tsolid");
  p->ioP0_.init(p->feSpaceP0_, "output_" + p->name_ + "/fieldsSolid");

  // assembly
  auto const & [assemblyIt, successEqn] =
      p->assemblies_.emplace(EQN_TYPE::CUSTOM, new AssemblySolid);
  assert(successEqn);
  p->equationType_ = EQN_TYPE::CUSTOM;

  auto assemblyPtr = dynamic_cast<AssemblySolid *>(assemblyIt->second.get());
  assemblyPtr->feSpaceGrad_.init(p->mesh_);
  assemblyPtr->feSpaceGradBC_.init(p->mesh_);
  assemblyPtr->uGrad_.init("uGradSolid", assemblyPtr->feSpaceGrad_);
  assemblyPtr->uGradBC_.init("uGradBCSolid", assemblyPtr->feSpaceGradBC_);
}

struct AssemblyFluid: public ProblemProXPDE::Assembly
{
  using Elem_T = ProblemProXPDEHeat::Elem_T;
  using FESpaceGradElem_T = proxpde::FESpace<
      ProblemProXPDEHeat::Mesh_T,
      ProblemProXPDEHeat::FE0_T::RefFE_T,
      ProblemProXPDEHeat::FE_T::RecommendedQR,
      Elem_T::dim>;
  static constexpr auto numPtsQRBC =
      proxpde::SideQR_T<typename ProblemProXPDEHeat::FESpace_T::QR_T>::numPts;
  using FESpaceGradBC_T = proxpde::FESpace<
      ProblemProXPDEHeat::Mesh_T,
      ProblemProXPDEHeat::FE0_T::RefFE_T,
      proxpde::SideGaussQR<ProblemProXPDEHeat::Elem_T, numPtsQRBC>,
      Elem_T::dim>;

  auto evaluate(ProblemProXPDE * pParent, proxpde::Builder<> & b) -> void override final
  {
    fmt::print("cht fluid assembly\n");
    auto p = dynamic_cast<ProblemProXPDEHeat *>(pParent);

    auto const alpha = p->params_["alpha"];

    // assembly lhs
    b.buildLhs(
        std::tuple{
            proxpde::AssemblyMass{1. / p->dt_, p->feSpace_},
            proxpde::AssemblyStiffness{alpha, p->feSpace_},
            proxpde::AssemblyAdvection{p->vel_, p->feSpace_},
        },
        p->bcs_);
    b.closeMatrix();

    // fmt::print("gradU min: {:.6e}\n", uGrad_.data.minCoeff());
    // fmt::print("gradU max: {:.6e}\n", uGrad_.data.maxCoeff());
    // assembly rhs
    b.buildRhs(
        std::tuple{
            proxpde::AssemblyProjection{1. / p->dt_, p->uOld_, p->feSpace_},
            // proxpde::AssemblyBCNaturalAnalytic{
            //     [](proxpde::Vec3 const &) { return proxpde::Vec1{1.5}; },
            //     proxpde::side::TOP,
            //     p->feSpace_},
            proxpde::AssemblyBCNaturalFE{
                uGradBC_, proxpde::Comp::y, proxpde::side::TOP, p->feSpace_},
        },
        p->bcs_);
  };

  FESpaceGradElem_T feSpaceGrad_;
  proxpde::FEVar<FESpaceGradElem_T> uGrad_;
  FESpaceGradBC_T feSpaceGradBC_;
  proxpde::FEVar<FESpaceGradBC_T> uGradBC_;
};

auto setupFluid(ProblemProXPDEHeat * p, proxpde::ParameterDict & config) -> void
{
  // name
  p->name_ = "proxpde_chts";

  // parameters
  auto const n = config["n"].as<uint>();
  config["mesh"]["type"] = proxpde::MeshType::STRUCTURED;
  config["mesh"]["origin"] = proxpde::Vec3{0.0, 0.0, 0.0};
  config["mesh"]["length"] = proxpde::Vec3{1.0, 0.75, 0.0};
  config["mesh"]["n"] = std::array{n, 3 * n / 4, 0u};
  config["mesh"]["flags"] = proxpde::Bitmask{proxpde::MeshFlags::NORMALS};

  // mesh
  proxpde::readMesh(p->mesh_, config["mesh"]);
  p->initMeshMED(p->name_, p->mesh_);

  // coupling
  p->couplingType_ = COUPLING_TYPE::MEDCOUPLING;

  // time
  p->time = 0.0;
  p->dt_ = config["dt"].as<double>();
  p->finalTime_ = config["final_time"].as<double>();
  p->printStep_ = config["print_step"].as<uint>();

  // finite elements
  p->feSpace_.init(p->mesh_);
  p->feSpaceP0_.init(p->mesh_);
  p->T_.init("T", p->feSpace_);
  p->T_ << 2.0;
  p->u_ = p->T_.data;

  // fields
  p->params_["alpha"] = 0.1;
  p->feSpaceVel_.init(p->mesh_);
  p->vel_.init("vel", p->feSpaceVel_);
  auto const velCoeff = config["vel_coeff"].as<double>();
  p->vel_ << [velCoeff](proxpde::Vec3 const & p) -> proxpde::Vec2
  {
    if (p[1] <= 0.75)
      return proxpde::Vec2{velCoeff * 4.0 * p[1] * (0.75 - p[1]), 0.0};
    return proxpde::Vec2{0.0, 0.0};
  };

  // bcs
  p->bcs_ = std::vector{
      proxpde::BCEss{
          p->feSpace_, proxpde::side::LEFT, [](proxpde::Vec3 const &) { return 1.0; }},
  };

  // io
  p->io_.init(p->feSpace_, "output_" + p->name_ + "/TFluid");
  p->ioP0_.init(p->feSpaceP0_, "output_" + p->name_ + "/fieldsFluid");

  p->initFieldMED("vel", "output_" + p->name_);
  p->setDataMED("vel", p->vel_.data, p->feSpaceVel_);

  // assembly
  auto const & [assemblyIt, successEqn] =
      p->assemblies_.emplace(EQN_TYPE::CUSTOM, new AssemblyFluid);
  assert(successEqn);
  p->equationType_ = EQN_TYPE::CUSTOM;

  auto assemblyPtr = dynamic_cast<AssemblyFluid *>(assemblyIt->second.get());
  assemblyPtr->feSpaceGrad_.init(p->mesh_);
  assemblyPtr->uGrad_.init("uGradFluid", assemblyPtr->feSpaceGrad_);
  assemblyPtr->uGrad_ << proxpde::Vec2{0.0, 1.5};
  assemblyPtr->feSpaceGradBC_.init(p->mesh_);
  assemblyPtr->uGradBC_.init("uGradFluid", assemblyPtr->feSpaceGradBC_);
  assemblyPtr->uGradBC_.data = assemblyPtr->uGrad_.data;
}

struct Mapper
{
  using Elem_T = ProblemProXPDEHeat::Elem_T;
  using Mesh_T = proxpde::Mesh<Elem_T>;
  using Mapping_T = std::unordered_map<proxpde::id_T, proxpde::id_T>;

  Mapper() = default;
  ~Mapper() = default;

  auto init(
      Mesh_T const & mesh1,
      proxpde::marker_T const marker1,
      Mesh_T const & mesh2,
      proxpde::marker_T const marker2) -> void
  {

    for (auto const & facet1: mesh1.facetList)
    {
      if (facet1.marker == marker1)
      {
        for (auto const & pt1: facet1.pts)
        {
          auto ptFound = false;
          for (auto const & facet2: mesh2.facetList)
          {
            if (facet2.marker == marker2)
            {
              for (auto const & pt2: facet2.pts)
              {
                if ((pt1->coord - pt2->coord).squaredNorm() < tol_ * tol_)
                {
                  if (ptMapping1to2_.contains(pt1->id))
                  {
                    assert(ptMapping1to2_.at(pt1->id) == pt2->id);
                    assert(ptMapping2to1_.at(pt2->id) == pt1->id);
                    // cannot check elemMapping, ptMapping can be set before
                    // assert(
                    //     elemMapping1to2_.at(facet1.facingElem[0].ptr->id) ==
                    //     facet2.facingElem[0].ptr->id);
                    // assert(
                    //     elemMapping2to1_.at(facet2.facingElem[0].ptr->id) ==
                    //     facet1.facingElem[0].ptr->id);
                  }
                  else
                  {
                    ptMapping1to2_[pt1->id] = pt2->id;
                    ptMapping2to1_[pt2->id] = pt1->id;
                    elemMapping1to2_[facet1.facingElem[0].ptr->id] =
                        facet2.facingElem[0].ptr->id;
                    elemMapping2to1_[facet2.facingElem[0].ptr->id] =
                        facet1.facingElem[0].ptr->id;
                  }
                  ptFound = true;
                }
                if (ptFound)
                  break;
              }
            }
            if (ptFound)
              break;
          }
        }
      }
    }
    // fmt::print("ptMapping1to2: {}\n", ptMapping1to2_);
    // fmt::print("elemMapping1to2: {}\n", elemMapping1to2_);
  }

  double const tol_ = 1.e-6;
  Mapping_T ptMapping1to2_;
  Mapping_T ptMapping2to1_;
  Mapping_T elemMapping1to2_;
  Mapping_T elemMapping2to1_;
};

struct Flux
{
  using Elem_T = ProblemProXPDEHeat::Elem_T;
  using Mesh_T = proxpde::Mesh<Elem_T>;
  using FE_T = ProblemProXPDEHeat::FE_T;
  using QR_T = FE_T::RecommendedQR;
  using QRSide_T =
      proxpde::SideGaussQR<proxpde::Quad, proxpde::SideQR<QR_T>::type::numPts>;
  using FESpaceSide_T =
      proxpde::FESpace<ProblemProXPDEHeat::Mesh_T, FE_T::RefFE_T, QRSide_T>;
  using FEFacet_T = FE_T::RefFE_T::FEFacet_T;
  using CurFEFacet_T =
      proxpde::CurFETraits<FEFacet_T, proxpde::SideQR<QR_T>::type>::type;

  Flux() = default;
  ~Flux() = default;

  auto init(
      proxpde::marker_T const marker,
      Mesh_T const & mesh,
      std::string_view name = "uSide",
      double const alpha = 1.0) -> void
  {
    assert(mesh.flags & proxpde::MeshFlags::NORMALS);
    marker_ = marker;
    feSpaceSide_.init(mesh);
    uSide_.init(name, feSpaceSide_);
    alpha_ = alpha;
  }

  auto compute(proxpde::Vec const & data) -> double
  {
    auto sum = 0.0;
    assert(data.size() == uSide_.data.size());
    uSide_.data = data;
    for (auto const & facet: feSpaceSide_.mesh->facetList)
    {
      if (facet.marker == marker_)
      {
        auto const [ePtr, side] = facet.facingElem[0];
        uSide_.reinit(*ePtr);
        curFEFacet_.reinit(facet);
        for (auto q = 0u; q < QRSide_T::numPtsSide; q++)
        {
          // auto const qpt = curFEFacet_.qpoint[q];
          // assert(std::fabs(qpt[1] - 0.75) < 1e-6);
          proxpde::Vec3 const gradSide =
              uSide_.evaluateGrad(QRSide_T::numPtsSide * side + q);
          sum += curFEFacet_.JxW[q] * alpha_ * gradSide.dot(facet._normal);
        }
      }
    }
    return sum;
  }

  proxpde::marker_T marker_ = proxpde::markerNotSet;
  FESpaceSide_T feSpaceSide_;
  proxpde::FEVar<FESpaceSide_T> uSide_;
  CurFEFacet_T curFEFacet_;
  double alpha_ = 1.0;
};

int main(int argc, char * argv[])
{
  proxpde::ParameterDict config;
  if (argc > 1)
  {
    config = YAML::LoadFile(argv[1]);
  }
  else
  {
    config["n"] = 20u;
    config["dt"] = 2.0e-3;
    config["final_time"] = 2.0;
    config["print_step"] = 100u;
    config["vel_coeff"] = 1.0;
  }
  config.validate({"n", "dt", "final_time", "print_step", "vel_coeff"});
  auto const n = config["n"].as<uint>();

  std::unique_ptr<ProblemProXPDEHeat> pSolid{new ProblemProXPDEHeat};
  setupSolid(pSolid.get(), config);
  pSolid->print();

  std::unique_ptr<ProblemProXPDEHeat> pFluid{new ProblemProXPDEHeat};
  setupFluid(pFluid.get(), config);
  pFluid->print();

  Flux fluxSolid;
  fluxSolid.init(
      proxpde::side::BOTTOM, pSolid->mesh_, "tSideSolid", pSolid->params_.at("alpha"));
  Flux fluxFluid;
  fluxFluid.init(
      proxpde::side::TOP, pFluid->mesh_, "tSideFluid", pFluid->params_.at("alpha"));

  Mapper mapper;
  mapper.init(pFluid->mesh_, proxpde::side::TOP, pSolid->mesh_, proxpde::side::BOTTOM);
  std::unordered_map<proxpde::DOFid_T, proxpde::id_T> ptMapISolid;
  for (auto k = 0u; k < pSolid->feSpace_.dof.ptMap.size(); k++)
    ptMapISolid[pSolid->feSpace_.dof.ptMap[k]] = k;
  std::unordered_map<proxpde::DOFid_T, proxpde::id_T> ptMapIFluid;
  for (auto k = 0u; k < pFluid->feSpace_.dof.ptMap.size(); k++)
    ptMapIFluid[pFluid->feSpace_.dof.ptMap[k]] = k;

  auto assemblySolid =
      dynamic_cast<AssemblySolid *>(pSolid->assemblies_.at(EQN_TYPE::CUSTOM).get());
  auto assemblyFluid =
      dynamic_cast<AssemblyFluid *>(pFluid->assemblies_.at(EQN_TYPE::CUSTOM).get());

  auto fluid2solid = std::unordered_map<proxpde::DOFid_T, proxpde::DOFid_T>{};
  for (auto dofTgt = 0u; dofTgt < assemblyFluid->uGrad_.data.size(); dofTgt++)
  {
    // dofTgt -> elemTgt -> elemSrc -> dofSrc
    // TODO: this relies on a specific DOF ordering
    auto const elemTgt = dofTgt / 2u;
    auto const compTgt = dofTgt % 2u;
    auto const & facetIds = pFluid->mesh_.elemToFacet[elemTgt];
    auto const markerFound =
        std::find_if(
            facetIds.begin(),
            facetIds.end(),
            [&pFluid](proxpde::id_T const & facetId)
            {
              if (facetId != proxpde::idNotSet)
                return pFluid->mesh_.facetList[facetId].marker == proxpde::side::TOP;
              return false;
            }) != facetIds.end();
    if (markerFound && compTgt == 1u)
    {
      auto const elemSrc = mapper.elemMapping1to2_[elemTgt];
      auto const dofSrc = 2 * elemSrc + compTgt;
      fluid2solid[dofTgt] = dofSrc;
    }
  }

  auto solid2fluid = std::unordered_map<proxpde::id_T, proxpde::DOFid_T>{};
  for (auto const & [dofTgt, posTgt]: pSolid->bcs_[0]._constrainedDofMap)
  {
    // dofTgt -> ptTgt -> ptSrc -> dofSrc
    auto const ptTgt = ptMapISolid.at(dofTgt);
    auto const ptSrc = mapper.ptMapping2to1_.at(ptTgt);
    auto const dofSrc = pFluid->feSpace_.dof.ptMap.at(ptSrc);
    solid2fluid[posTgt] = dofSrc;
  }

  auto const h = 1.0 / n;

  auto const externalIters = 1u;
  auto const internalItersSolid = 1u;
  auto const internalItersFluid = 1u;
  // bool updateDt = true;
  while (pSolid->run() || pFluid->run())
  {
    // if (updateDt && pFluid->time > 1.0 - 1e-6)
    // {
    //   updateDt = false;
    //   pFluid->dt_ *= 2;
    //   pFluid->printStep_ /= 2u;
    //   pSolid->dt_ *= 2;
    //   pSolid->printStep_ /= 2u;
    // }

    pFluid->advance();
    pSolid->advance();

    for (uint ei = 0u; ei < externalIters; ei++)
    {
      // update bc: flux solid (src) -> fluid (tgt)
      proxpde::computeGradient(
          assemblySolid->uGrad_.data,
          *assemblySolid->uGrad_.feSpace,
          pSolid->T_.data,
          *pSolid->T_.feSpace);
      auto flux = 0.0;
      for (auto const & [dofTgt, dofSrc]: fluid2solid)
      {
        auto const fluxLocal = assemblySolid->uGrad_.data[dofSrc];
        assemblyFluid->uGrad_.data[dofTgt] = fluxLocal;
        flux += fluxLocal * h;
      }
      assemblyFluid->uGradBC_.data = assemblyFluid->uGrad_.data;
      // assemblyFluid->uGradBC_ << proxpde::Vec2{0.0, 1.5};

      fmt::print("total flux s->f: {:.6e}\n", flux);
      if (std::fabs(flux) > 10.)
        std::abort();

      // solve fluid
      for (auto ii = 0u; ii < internalItersFluid; ii++)
        pFluid->solve();

      // update bc: value fluid (src) -> solid (tgt)
      for (auto const & [posTgt, dofSrc]: solid2fluid)
        pSolid->bcs_[0].data[posTgt] = pFluid->T_.data[dofSrc];

      // solve solid
      for (auto ii = 0u; ii < internalItersSolid; ii++)
        pSolid->solve();

      // compute flux from solid
      auto const currentFluxSolid = fluxSolid.compute(pSolid->T_.data);
      fmt::print("fluxSolid: {:.6e}\n", currentFluxSolid);

      // compute flux to fluid
      auto const currentFluxFluid = fluxFluid.compute(pFluid->T_.data);
      fmt::print("fluxFluid: {:.6e}\n", currentFluxFluid);
    }

    pFluid->print();
    pSolid->print();
  }

  return 0;
}
