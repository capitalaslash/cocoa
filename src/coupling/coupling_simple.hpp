#pragma once

// std
#include <vector>

// local
#include "coupling/coupling_manager.hpp"
#include "coupling/mesh_coupling.hpp"

namespace cocoa
{

struct CouplingSimple: public CouplingManager
{
  struct Matrix
  {
    Matrix() = default;
    Matrix(size_t n, size_t m): data(n, std::vector<double>(m)) {}
    ~Matrix() = default;

    void init(size_t n, size_t m) { data.resize(n, std::vector<double>(m)); }

    double & operator()(size_t i, size_t j) { return data[i][j]; }

    std::vector<std::vector<double>> data;
  };

  CouplingSimple(COUPLING_SCOPE scope): CouplingManager(COUPLING_TYPE::SIMPLE, scope)
  {
    if (scope == COUPLING_SCOPE::BOUNDARY)
    {
      fmt::println(stderr, "boundary coupling not yet implemented!");
      std::abort();
    }
  }
  ~CouplingSimple() = default;

  void setup(CouplingInterface interfaceSrc, CouplingInterface interfaceTgt) override;
  void initFieldCoupling() override;

  void project(std::string_view fieldNameSrc, std::string_view fieldNameTgt) override;

  Matrix m_;
};

} // namespace cocoa
