#pragma once

#include "plugins.hpp"

#ifdef COCOA_ENABLE_PROXPDE

// std
#include <vector>

// proxpde
#include <proxpde/def.hpp>

#include <proxpde/bc.hpp>
#include <proxpde/fe.hpp>
#include <proxpde/fespace.hpp>
#include <proxpde/geo.hpp>
#include <proxpde/iomanager.hpp>
#include <proxpde/mesh.hpp>
#include <proxpde/var.hpp>

// local
#include "enums.hpp"
#include "problem/problem.hpp"

namespace cocoa
{

struct ProblemProXPDE: public Problem
{
  struct Assembly
  {
    Assembly() = default;
    virtual ~Assembly() = default;
    virtual auto evaluate(ProblemProXPDE *, proxpde::Builder<> &) -> void = 0;
  };

  using Elem_T = proxpde::Quad;
  using Mesh_T = proxpde::Mesh<Elem_T>;

  ProblemProXPDE();
  virtual ~ProblemProXPDE() = default;

  virtual void setup(Problem::ConfigList_T const & configs) override = 0;
  void advance() override final;
  bool run() const override final;
  virtual uint solve() override;
  virtual void print() override = 0;

  virtual std::vector<std::string> varNames() override = 0;
  Marker findRegion(std::string_view name) override;

  virtual uint size() const = 0;

  std::unique_ptr<MeshCoupling> initMeshCoupling(COUPLING_TYPE type) override;
  virtual std::unique_ptr<FieldCoupling> initFieldCoupling(
      COUPLING_TYPE type,
      std::string_view name,
      MeshCoupling const * mesh) override = 0;
  virtual void setFieldData(FieldCoupling * field) override = 0;
  virtual void getFieldData(FieldCoupling const & field) override = 0;
  // void initFieldMED(std::string_view fieldName, std::string_view path);

  // template <typename FESpace>
  // void setDataMED(
  //     std::string_view fieldName, proxpde::Vec const & u, FESpace const & feSpace);

  // template <typename FESpace>
  // void
  // getDataMED(std::string_view fieldName, proxpde::Vec & u, FESpace const & feSpace);

  std::string name_;
  EQN_TYPE equationType_;
  Mesh_T mesh_;
  std::unordered_map<EQN_TYPE, std::unique_ptr<Assembly>> assemblies_;

  std::vector<std::string> couplingExport_;
  std::vector<std::string> couplingImport_;
  proxpde::Vec u_;
  proxpde::Vec uOld_;
  double dt_;
  double finalTime_;
  std::filesystem::path outputPrefix_ = "./tmp";

  static std::unique_ptr<Problem> build(EQN_TYPE const type);
};

// =====================================================================

struct ProblemProXPDEHeat: public ProblemProXPDE
{
  using FE_T = proxpde::LagrangeFE<Elem_T, 1u>;
  using FESpace_T = proxpde::FESpace<Mesh_T, FE_T::RefFE_T, FE_T::RecommendedQR>;
  using FE0_T = proxpde::LagrangeFE<Elem_T, 0u>;
  using FESpaceP0_T = proxpde::FESpace<Mesh_T, FE0_T::RefFE_T, FE_T::RecommendedQR>;
  using FESpaceVel_T = proxpde::FESpace<Mesh_T, FE_T::RefFE_T, FE_T::RecommendedQR, 2u>;

  ProblemProXPDEHeat();
  ~ProblemProXPDEHeat() = default;

  void setup(Problem::ConfigList_T const & configs) override;
  uint solve() override;
  void print() override;

  std::unique_ptr<FieldCoupling> initFieldCoupling(
      COUPLING_TYPE type, std::string_view name, MeshCoupling const * mesh) override;
  void setFieldData(FieldCoupling * field) override;
  void getFieldData(FieldCoupling const & field) override;

  std::vector<std::string> varNames() override { return {"T"}; }

  virtual uint size() const override { return feSpace_.dof.size; }

  FESpace_T feSpace_;
  proxpde::FEVar<FESpace_T> T_;
  FESpaceP0_T feSpaceP0_;
  std::unordered_map<std::string, double> params_;
  std::unordered_map<std::string, proxpde::FEVar<FESpace_T>> fieldsP1_;
  std::unordered_map<std::string, proxpde::FEVar<FESpaceP0_T>> fieldsP0_;
  FESpaceVel_T feSpaceVel_;
  proxpde::FEVar<FESpaceVel_T> vel_;
  std::vector<proxpde::BCEss<FESpace_T>> bcs_;
  proxpde::IOManager<FESpace_T> io_;
  proxpde::IOManager<FESpaceP0_T> ioP0_;
};

// =====================================================================

struct ProblemProXPDENS: public ProblemProXPDE
{
  using Elem_T = proxpde::Quad;
  using Mesh_T = proxpde::Mesh<Elem_T>;
  using FE0_T = proxpde::LagrangeFE<Elem_T, 0u>;
  using FE1_T = proxpde::LagrangeFE<Elem_T, 1u>;
  using FE2_T = proxpde::LagrangeFE<Elem_T, 2u>;
  using FESpaceP_T = proxpde::FESpace<Mesh_T, FE1_T::RefFE_T, FE2_T::RecommendedQR>;
  using FESpaceVel_T =
      proxpde::FESpace<Mesh_T, FE2_T::RefFE_T, FE2_T::RecommendedQR, 2u>;
  using FESpaceVelQ1_T =
      proxpde::FESpace<Mesh_T, FE1_T::RefFE_T, FE2_T::RecommendedQR, 2u>;
  using FESpaceP0_T = proxpde::FESpace<Mesh_T, FE0_T::RefFE_T, FE2_T::RecommendedQR>;

  ProblemProXPDENS();
  ~ProblemProXPDENS() = default;

  void setup(Problem::ConfigList_T const & configs) override;
  uint solve() override;
  void print() override;

  std::unique_ptr<FieldCoupling> initFieldCoupling(
      COUPLING_TYPE type, std::string_view name, MeshCoupling const * mesh) override;
  void setFieldData(FieldCoupling * field) override;
  void getFieldData(FieldCoupling const & field) override;

  std::vector<std::string> varNames() override { return {"p", "vel"}; }

  virtual uint size() const override
  {
    return FESpaceVel_T::dim * feSpaceVel_.dof.size + feSpaceP_.dof.size;
  }

  FESpaceVel_T feSpaceVel_;
  proxpde::FEVar<FESpaceVel_T> vel_;
  FESpaceP_T feSpaceP_;
  proxpde::FEVar<FESpaceP_T> p_;
  std::unordered_map<std::string, proxpde::FEVar<FESpaceP_T>> fieldsP1_;
  FESpaceP0_T feSpaceP0_;
  std::unordered_map<std::string, proxpde::FEVar<FESpaceP0_T>> fieldsP0_;
  FESpaceVelQ1_T feSpaceVelQ1_;
  proxpde::FEVar<FESpaceVelQ1_T> velQ1_;
  proxpde::L2Projector<FESpaceVelQ1_T, FESpaceVel_T> projectorQ2Q1_;
  double viscosity_;
  std::vector<proxpde::BCEss<FESpaceP_T>> bcsP_;
  std::vector<proxpde::BCEss<FESpaceVel_T>> bcsVel_;
  proxpde::IOManager<FESpaceVel_T> ioVel_;
  proxpde::IOManager<FESpaceP_T> ioP_;
};

} // namespace cocoa

#endif
