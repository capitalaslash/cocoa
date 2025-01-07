#pragma once

#include "plugins.hpp"

#ifdef COCOA_ENABLE_MEDCOUPLING

// std
#include <string_view>

// medcoupling
#include <MEDCouplingFieldDouble.hxx>

// local
#include "coupling/field_coupling.hpp"
#include "enums.hpp"

inline MEDCoupling::TypeOfField supportType2MEDtype(SUPPORT_TYPE const support)
{
  switch (support)
  {
  case SUPPORT_TYPE::ON_NODES:
    return MEDCoupling::ON_NODES;
  case SUPPORT_TYPE::ON_CELLS:
    return MEDCoupling::ON_CELLS;
  default:
    fmt::print(stderr, "field type not recognized\n");
    return MEDCoupling::ON_NODES; // no default option available in MEDCoupling
  }
}

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
  double const * dataPtr() const override
  {
    return fieldPtr_->getArray()->getConstPointer();
  }
  double operator[](size_t k) const override { return *(this->dataPtr() + k); }

  // std::vector<double> getData() override;
  void
  init(std::string_view name, MeshCoupling * mesh, SUPPORT_TYPE const support) override;
  void setValues(std::span<double> const & data, uint dim = 1U) override;
  void setValues(double value, uint size, uint dim = 1U) override;

  void printVTK(double time, uint iter) override;
  void printPVD() const;

  MEDCoupling::MEDCouplingFieldDouble * fieldPtr_;
  bool inited_ = false;
  std::vector<Frame> frames;
};

#endif
