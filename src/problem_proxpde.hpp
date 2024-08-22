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
#include "field_coupling.hpp"
#include "problem.hpp"

struct ProblemProXPDE: public Problem
{
  using ParamList_T = Problem::ParamList_T;
  using Elem_T = proxpde::Quad;
  using Mesh_T = proxpde::Mesh<Elem_T>;
  using FE_T = proxpde::LagrangeFE<Elem_T, 1U>;
  using FESpace_T = proxpde::FESpace<Mesh_T, FE_T::RefFE_T, FE_T::RecommendedQR>;

  ProblemProXPDE(): Problem{PROBLEM_TYPE::PROXPDE} {}
  ~ProblemProXPDE() = default;

  void setup(ParamList_T const & params) override;
  void initMeshMED(std::string_view meshName);
  void initFieldMED(std::string_view fieldName);
  void updateFieldMED(std::string_view fieldName);
  void advance() override;
  bool run() override;
  void solve() override;
  void print() override;

  Mesh_T mesh_;
  FESpace_T feSpace_;
  std::vector<proxpde::BCEss<FESpace_T>> bcs_;
  proxpde::IOManager<FESpace_T> io_;
  double diff_;
  proxpde::FEVar<FESpace_T> u_;
  proxpde::Vec uOld_;
  double dt_;
  double finalTime_;
};

#endif
