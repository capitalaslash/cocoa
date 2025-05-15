// local
#include "problem/problem_proxpde.hpp"

static proxpde::marker_T constexpr LEFT_BOTTOM = 101u;
static proxpde::marker_T constexpr LEFT_TOP = 102u;
static proxpde::marker_T constexpr INTERFACE = 201u;

using namespace cocoa;

struct AssemblyCHTM: public ProblemProXPDE::Assembly
{
  auto evaluate(ProblemProXPDE * pParent, proxpde::Builder<> & b) -> void override
  {
    fmt::println("proxpde chtm assembly");
    auto p = dynamic_cast<ProblemProXPDEHeat *>(pParent);

    auto & alpha = p->fieldsP0_.at("alpha");
    auto & q = p->fieldsP0_.at("q");

    // assembly lhs
    b.buildLhs(
        std::tuple{
            proxpde::AssemblyMass{1. / p->dt_, p->feSpace_},
            proxpde::AssemblyStiffnessFE{alpha, p->feSpace_},
            // proxpde::AssemblyStiffness{p->alpha_, p->feSpace_},
            proxpde::AssemblyAdvection{p->vel_, p->feSpace_},
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
};

auto setup(ProblemProXPDEHeat * p, proxpde::ParameterDict & config) -> void
{
  // name
  p->name_ = "proxpde_chtm";

  // parameters
  auto const n = config["n"].as<uint>();
  config["mesh"]["type"] = proxpde::MeshType::STRUCTURED;
  config["mesh"]["origin"] = proxpde::Vec3{0.0, 0.0, 0.0};
  config["mesh"]["length"] = proxpde::Vec3{1.0, 1.0, 0.0};
  config["mesh"]["n"] = std::array{n, n, 0u};
  config["mesh"]["flags"] =
      proxpde::MeshFlags::INTERNAL_FACETS | proxpde::MeshFlags::NORMALS;

  // mesh
  proxpde::readMesh(p->mesh_, config["mesh"]);

  // auto const k = 0.25;
  // auto const a = std::pow(k, 1.0 / (N - 1));
  // // auto const hGeom = 1.0 * (1.0 - a) / (1 - k * a);
  // auto const hFixed = 1.0 / N;
  // for (auto & pt: p->mesh_.pointList)
  // {
  //   uint const i = std::round(pt.coord[0] / hFixed);
  //   uint const j = std::round(pt.coord[1] / hFixed);
  //   pt.coord = proxpde::Vec3{
  //       1.0 - 1.0 * (1.0 - std::pow(a, N - i)) / (1.0 - k * a),
  //       1.0 * (1.0 - std::pow(a, j)) / (1.0 - k * a),
  //       0.0};
  // }

  // auto const hFixed = 1.0 / N;
  // for (auto & pt: p->mesh_.pointList)
  // {
  //   uint const i = std::round(pt.coord[0] / hFixed);
  //   uint const j = std::round(pt.coord[1] / hFixed);

  //   auto dx = i * hFixed;
  //   if (i < N / 2)
  //     dx = 0.0 + i * 0.5 * hFixed;
  //   else
  //     dx = 0.25 + (i - N / 2) * 1.5 * hFixed;

  //   auto dy = j * hFixed;
  //   if (j < 5 * N / 12)
  //     dy = 0.0 + j * 1.5 * hFixed;
  //   else if (j < 11 * N / 12)
  //     dy = 0.625 + (j - 5 * N / 12) * 0.5 * hFixed;
  //   else
  //     dy = 0.875 + (j - 11 * N / 12) * 1.5 * hFixed;

  //   pt.coord = proxpde::Vec3{dx, dy, 0.0};
  // }

  for (auto & facet: p->mesh_.facetList)
  {
    if (std::fabs(facet.midpoint()[1] - 0.75) < 1.0e-6)
    {
      facet.marker = INTERFACE;
    }
  }

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
  p->fieldsP0_.emplace("alpha", proxpde::FEVar{"alpha", p->feSpaceP0_});
  p->fieldsP0_.at("alpha") << [](proxpde::Vec3 const & p) -> double
  {
    if (p[1] > 0.75)
      return 1.0;
    return 0.1;
  };
  p->fieldsP0_.emplace("q", proxpde::FEVar{"q", p->feSpaceP0_});
  p->fieldsP0_.at("q") << [](proxpde::Vec3 const & p) -> double
  {
    if (p[1] > 0.75)
      return 6.0;
    return 0.0;
  };
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
  for (auto & facet: p->mesh_.facetList)
  {
    if (facet.marker == proxpde::side::LEFT)
    {
      if (facet.midpoint()[1] < 0.75)
        facet.marker = LEFT_BOTTOM;
      else
        facet.marker = LEFT_TOP;
    }
  }
  auto const one = [](proxpde::Vec3 const &) { return 1.0; };
  p->bcs_ = std::vector{
      proxpde::BCEss{p->feSpace_, LEFT_BOTTOM, one},
  };

  // io
  p->io_.init(p->feSpace_, "output_" + p->name_ + "/T");
  p->ioP0_.init(p->feSpaceP0_, "output_" + p->name_ + "/fields");

  p->initFieldMED("vel", "output_" + p->name_);
  p->setDataMED("vel", p->vel_.data, p->feSpaceVel_);

  auto const & [_, successEqn] =
      p->assemblies_.emplace(EQN_TYPE::CUSTOM, new AssemblyCHTM);
  assert(successEqn);
  p->equationType_ = EQN_TYPE::CUSTOM;
}

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
    config["dt"] = 1.0e+0;
    config["final_time"] = 10.0;
    config["print_step"] = 2u;
    config["vel_coeff"] = 1.0;
  }
  config.validate({"n", "dt", "final_time", "print_step", "vel_coeff"});

  std::unique_ptr<ProblemProXPDEHeat> p{new ProblemProXPDEHeat};

  // p->setup({{"config_file", "proxpde_heat.yaml"}});
  setup(p.get(), config);

  p->print();

  using FE_T = ProblemProXPDEHeat::FE_T;
  using QR_T = FE_T::RecommendedQR;
  using QRSide_T =
      proxpde::SideGaussQR<proxpde::Quad, proxpde::SideQR<QR_T>::type::numPts>;
  using FESpaceSide_T =
      proxpde::FESpace<ProblemProXPDEHeat::Mesh_T, FE_T::RefFE_T, QRSide_T>;
  using FEFacet_T = FE_T::RefFE_T::FEFacet_T;
  using CurFEFacet_T =
      proxpde::CurFETraits<FEFacet_T, proxpde::SideQR<QR_T>::type>::type;
  FESpaceSide_T feSpaceGrad{p->mesh_};
  proxpde::FEVar tSide{"Tside", feSpaceGrad};
  auto & alpha = p->fieldsP0_.at("alpha");
  CurFEFacet_T curFEFacet;

  while (p->run())
  {
    p->advance();
    p->solve();
    p->print();

    // compute flux across interface
    tSide.data = p->T_.data;
    std::unordered_set<proxpde::id_T> processedFacets;
    auto fluxSolid = 0.0;
    auto fluxFluid = 0.0;
    for (auto const & facet: p->mesh_.facetList)
    {
      if (facet.marker == INTERFACE)
      {
        for (auto const & [ePtr, side]: facet.facingElem)
        {
          tSide.reinit(*ePtr);
          alpha.reinit(*ePtr);
          curFEFacet.reinit(facet);
          double const alphaLocal = alpha.dataLocal[0];

          // solid side
          if (ePtr->midpoint()[1] > 0.75)
          {
            for (auto q = 0u; q < QRSide_T::numPtsSide; q++)
            {
              auto const qpt = curFEFacet.qpoint[q];
              assert(std::fabs(qpt[1] - 0.75) < 1e-6);
              proxpde::Vec3 const gradSide =
                  tSide.evaluateGrad(QRSide_T::numPtsSide * side + q);
              fluxSolid += curFEFacet.JxW[q] * alphaLocal * gradSide[1];
            }
          }
          // fluid side
          else
          {
            for (auto q = 0u; q < QRSide_T::numPtsSide; q++)
            {
              auto const qpt = curFEFacet.qpoint[q];
              assert(std::fabs(qpt[1] - 0.75) < 1e-6);
              proxpde::Vec3 const gradSide =
                  tSide.evaluateGrad(QRSide_T::numPtsSide * side + q);
              fluxFluid += curFEFacet.JxW[q] * alphaLocal * gradSide[1];
            }
          }
        }
      }
    }
    fmt::println("fluxSolid: {:.6e}", fluxSolid);
    fmt::println("fluxFluid: {:.6e}", fluxFluid);
  }

  return 0;
}
