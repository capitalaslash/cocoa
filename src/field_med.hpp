#pragma once

// std
#include <string_view>

// medcoupling
#include "MEDCouplingFieldDouble.hxx"

// local
#include "field_coupling.hpp"
#include "med_mesh.hpp"

struct FieldMED: public FieldCoupling
{
  struct Frame
  {
    uint iter;
    double time;
  };

  FieldMED() = default;
  ~FieldMED();

  double * getData() override { return nullptr; }
  void init(std::string_view name, MEDMesh & mesh) override;
  void initIO(std::string_view filename) override;
  void setValues(std::vector<double> const & data) override;
  void setValues(double value, uint size) override;

  void printVTK(double time, uint iter);
  void printPVD();

  MEDCoupling::MEDCouplingFieldDouble * fieldPtr_;
  bool inited_ = false;
  std::vector<Frame> frames;
};
