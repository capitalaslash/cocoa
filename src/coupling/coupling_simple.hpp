#pragma once

// std
#include <vector>

// local
#include "coupling/coupling_manager.hpp"

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

  CouplingSimple(): CouplingManager(COUPLING_TYPE::SIMPLE) {}
  ~CouplingSimple() = default;

  void setup(Problem * pSrc, Problem * pTgt) override;

  void project(std::string_view fieldNameSrc, std::string_view fieldNameTgt) override;

  Matrix m_;
};
