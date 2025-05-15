// local
#include "problem/problem_proxpde.hpp"

auto constexpr N = 20u;

using namespace cocoa;

struct AssemblyCustom: public ProblemProXPDE::Assembly
{
  auto evaluate(ProblemProXPDE * pParent, proxpde::Builder<> & b) -> void override
  {
    fmt::println("proxpde custom assembly");
    auto p = dynamic_cast<ProblemProXPDEHeat *>(pParent);
    auto const alpha = p->params_.at("alpha");
    auto const & q = p->fieldsP0_.at("q");

    // assembly lhs
    b.buildLhs(
        std::tuple{
            proxpde::AssemblyMass{1. / p->dt_, p->feSpace_},
            proxpde::AssemblyStiffness{alpha, p->feSpace_},
            proxpde::AssemblyAdvection{p->vel_, p->feSpace_},
        },
        p->bcs_);
    b.closeMatrix();

    // assembly rhs
    b.buildRhs(
        std::tuple{
            proxpde::AssemblyProjection{1. / p->dt_, p->uOld_, p->feSpace_},
            proxpde::AssemblyProjection{1., q.data, p->feSpaceP0_, p->feSpace_},
            proxpde::AssemblyBCNaturalAnalytic{
                [](proxpde::Vec3 const & p) { return proxpde::Vec1{1.5}; },
                proxpde::side::TOP,
                p->feSpace_}},
        p->bcs_);
  }
};

auto setup(ProblemProXPDEHeat * p) -> void
{
  // name
  p->name_ = "proxpde_custom";

  // parameters
  proxpde::ParameterDict config;
  config["mesh"]["type"] = proxpde::MeshType::STRUCTURED;
  config["mesh"]["origin"] = proxpde::Vec3{0.0, 0.0, 0.0};
  config["mesh"]["length"] = proxpde::Vec3{1.0, 0.75, 0.0};
  config["mesh"]["n"] = std::array{N, 3 * N / 4, 0u};

  // mesh
  proxpde::readMesh(p->mesh_, config["mesh"]);
  p->initMeshMED(p->name_, p->mesh_);

  // coupling
  p->couplingType_ = COUPLING_TYPE::MEDCOUPLING;

  // time
  p->time = 0.0;
  p->dt_ = 1.0;
  p->finalTime_ = 10.0;

  // finite elements
  p->feSpace_.init(p->mesh_);
  p->feSpaceP0_.init(p->mesh_);
  p->T_.init("T", p->feSpace_);
  p->T_ << 0.0;
  p->u_ = p->T_.data;

  // fields
  p->params_["alpha"] = 0.1;
  p->fieldsP0_.emplace("q", proxpde::FEVar{"q", p->feSpaceP0_});
  p->fieldsP0_.at("q") << 0.0;
  p->feSpaceVel_.init(p->mesh_);
  p->vel_.init("vel", p->feSpaceVel_);
  p->vel_ << [](proxpde::Vec3 const & p)
  { return proxpde::Vec2{4.0 * p[1] * (0.75 - p[1]), 0.0}; };

  // bcs
  p->bcs_ = std::vector{
      proxpde::BCEss{
          p->feSpace_, proxpde::side::LEFT, [](proxpde::Vec3 const &) { return 1.0; }},
  };

  // io
  p->io_.init(p->feSpace_, "output_" + p->name_ + "/T");
  p->ioP0_.init(p->feSpaceP0_, "output_" + p->name_ + "/fields");

  // coupling - again
  p->couplingExport_.push_back("Tcfd");
  p->initFieldMED(p->couplingExport_[0], "output_" + p->name_);
  p->setDataMED(p->couplingExport_[0], p->T_.data, p->feSpace_);

  p->initFieldMED("vel", "output_" + p->name_);
  p->setDataMED("vel", p->vel_.data, p->feSpaceVel_);

  // assembly
  [[maybe_unused]] auto const & [_, success] =
      p->assemblies_.emplace(EQN_TYPE::CUSTOM, new AssemblyCustom);
  assert(success);
  p->equationType_ = EQN_TYPE::CUSTOM;
}

int main()
{
  std::unique_ptr<ProblemProXPDEHeat> p{new ProblemProXPDEHeat};

  // p->setup({{"config_file", "proxpde_heat.yaml"}});
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
