#pragma once

// std
#include <cassert>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// local
#include "field_coupling.hpp"
#include "mesh_coupling.hpp"
#include "problem.hpp"

void setAssemblies();

struct ProblemFD1D: public Problem
{
  using Assembly_T = std::function<void(ProblemFD1D *)>;

  struct Matrix
  {
    Matrix() = default;
    explicit Matrix(uint n): diag(n), diagUp(n), diagDown(n) {}
    ~Matrix() = default;

    void init(uint n)
    {
      diag.resize(n);
      diagUp.resize(n);
      diagDown.resize(n);
    }

    std::vector<double> diag;
    std::vector<double> diagUp;
    std::vector<double> diagDown;
  };

  ProblemFD1D() = default;
  ~ProblemFD1D() = default;

  void setup(Problem::ParamList_T const & params) override;
  bool run() override;
  void advance() override;
  void solve() override;
  void print() override;

  void initMeshMED(std::filesystem::path const & fileName);
  void initFieldMED(std::filesystem::path const & fileName);

  FieldCoupling getField(std::string_view name) override { return uCoupling_; }
  void setField(std::string_view name, FieldCoupling const & field) override
  {
    assert(name == uCoupling_.name_);
    uCoupling_ = field;
  }

  double start_;
  double h_;
  uint n_;
  MeshCoupling meshCoupling_;
  std::vector<double> u_;
  FieldCoupling uCoupling_;
  std::vector<double> uOld_;
  std::vector<double> q_;
  double diff_;
  double finalTime_;
  double dt_;
  Matrix m_;
  std::vector<double> rhs_;
  std::string assemblyName_;
  double bcStart_;
  double bcEnd_;
  std::string outFile_ = "fd1d";
  std::string nameExt_ = "uExternal";

  static std::unordered_map<std::string, Assembly_T> assemblies_;
};
