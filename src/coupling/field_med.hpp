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

namespace cocoa
{

inline constexpr MEDCoupling::TypeOfField supportType2MED(SUPPORT_TYPE const support)
{
  using enum SUPPORT_TYPE;
  switch (support)
  {
  case ON_NODES:
    return MEDCoupling::ON_NODES;
  case ON_CELLS:
    return MEDCoupling::ON_CELLS;
  default:
    fmt::println(stderr, "field support not recognized");
    // // requires c++23
    // fmt::println(
    //     stderr, "field support {} not recognized", std::to_underlying(support));
    std::abort();
    return MEDCoupling::ON_NODES; // no default option available in MEDCoupling
  }
}

inline constexpr MEDCoupling::NatureOfField natureType2MED(NATURE_TYPE const nature)
{
  using enum NATURE_TYPE;
  switch (nature)
  {
  case INTENSIVE_MAXIMUM:
    return MEDCoupling::IntensiveMaximum;
  case EXTENSIVE_MAXIMUM:
    return MEDCoupling::ExtensiveMaximum;
  case INTENSIVE_CONSERVATION:
    return MEDCoupling::IntensiveConservation;
  case EXTENSIVE_CONSERVATION:
    return MEDCoupling::ExtensiveConservation;
  default:
    fmt::println(stderr, "field nature not recognized");
    std::abort();
    return MEDCoupling::NoNature;
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
  double at(size_t k) const override { return *(this->dataPtr() + k); }
  double operator[](size_t k) const override { return *(this->dataPtr() + k); }

  // std::vector<double> getData() override;
  void init(
      std::string_view name,
      MeshCoupling const * mesh,
      SUPPORT_TYPE support,
      NATURE_TYPE nature) override;
  void setValues(std::span<const double> data, uint dim = 1u) override;
  void setValues(double value, uint size, uint dim = 1u) override;

  void printVTK(double time, uint iter) override;
  void printPVD() const;

  MEDCoupling::MEDCouplingFieldDouble * fieldPtr_;
  bool inited_ = false;
  std::vector<Frame> frames;
};

} // namespace cocoa

#endif
