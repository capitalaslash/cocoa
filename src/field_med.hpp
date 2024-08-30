#pragma once

// std
#include <string_view>

// medcoupling
#include "MEDCouplingFieldDouble.hxx"

// local
#include "enums.hpp"
#include "field_coupling.hpp"

struct FieldMED: public FieldCoupling
{
  struct Frame
  {
    uint iter;
    double time;
  };

  FieldMED(): FieldCoupling{COUPLING_TYPE::MEDCOUPLING} {}
  ~FieldMED();

  size_t size() const noexcept override { return fieldPtr_->getNumberOfValues(); }
  double const * dataPtr() override { return fieldPtr_->getArray()->getConstPointer(); }
  // std::vector<double> getData() override;
  void init(std::string_view name, MeshCoupling * mesh) override;
  void initIO(std::string_view filename) override;
  void setValues(std::vector<double> const & data, uint dim = 1U) override;
  void setValues(double value, uint size, uint dim = 1U) override;

  void printVTK(double time, uint iter) override;
  void printPVD() const;

  MEDCoupling::MEDCouplingFieldDouble * fieldPtr_;
  bool inited_ = false;
  std::vector<Frame> frames;
};
