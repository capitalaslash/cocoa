#pragma once

// std
#include <filesystem>
#include <memory>
#include <span>
#include <string>

// fmt
#include <fmt/core.h>

// local
#include "coupling/mesh_coupling.hpp"
#include "enums.hpp"

enum struct SUPPORT_TYPE : uint8_t
{
  NONE = 0,
  ON_NODES,
  ON_CELLS,
};

struct FieldCoupling
{
  FieldCoupling() = default;
  explicit FieldCoupling(COUPLING_TYPE type): type_{type} {}
  virtual ~FieldCoupling() = default;

  virtual size_t size() const noexcept = 0;
  virtual double const * dataPtr() const = 0;
  virtual double operator[](size_t k) const = 0;
  virtual double at(size_t k) const = 0;
  // virtual std::vector<double> getData() = 0;
  virtual void
  init(std::string_view name, MeshCoupling * mesh, SUPPORT_TYPE const support) = 0;
  virtual void initIO(std::filesystem::path const & prefix)
  {
    prefix_ = prefix;
    std::filesystem::create_directories(prefix_);
  }
  virtual void setValues(std::span<double> const & data, uint const dim = 1U) = 0;
  virtual void setValues(double value, uint size, uint const dim = 1U) = 0;
  virtual void printVTK(double time, uint iter) = 0;

  static std::unique_ptr<FieldCoupling> build(COUPLING_TYPE type);
  static std::unique_ptr<FieldCoupling> build(std::string_view type);

  std::string name_ = "";
  COUPLING_TYPE type_ = COUPLING_TYPE::NONE;
  std::filesystem::path prefix_ = "./tmp";
};

inline std::ostream & operator<<(std::ostream & out, FieldCoupling const & f)
{
  out << "FieldCoupling(" << f.name_ << ", " << static_cast<uint>(f.type_) << ")";
  return out;
}

// =====================================================================

struct FieldSimple: public FieldCoupling
{
  FieldSimple(): FieldCoupling{COUPLING_TYPE::SIMPLE} {}
  ~FieldSimple() = default;

  size_t size() const noexcept override { return data_.size(); }
  double const * dataPtr() const override { return data_.data(); }
  double at(size_t k) const override { return data_[k]; }
  double operator[](size_t k) const override { return data_[k]; }

  // std::span<double> getData() override { return data_; }
  virtual void init(
      std::string_view name,
      MeshCoupling * /*mesh*/,
      SUPPORT_TYPE const supportType) override
  {
    name_ = name;
    supportType_ = supportType;
  }
  virtual void initIO(std::filesystem::path const & /*prefix*/) override {}
  virtual void
  setValues(std::span<double> const & data, uint const /*dim*/ = 1U) override
  {
    data_.resize(data.size());
    std::copy(data.begin(), data.end(), data_.begin());
  }
  virtual void setValues(double value, uint size, uint const /*dim*/ = 1U) override
  {
    data_.resize(size, value);
  }

  virtual void printVTK(double /*time*/, uint /*iter*/) override
  {
    if (messageVTK_)
    {
      fmt::print(stderr, "VTK print is not implemented in FieldSimple.\n");
      messageVTK_ = false;
    }
  }

  std::vector<double> data_;
  SUPPORT_TYPE supportType_;
  bool messageVTK_ = true;
};
