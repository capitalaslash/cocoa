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

namespace cocoa
{

enum struct SUPPORT_TYPE : uint8_t
{
  NONE = 0u,
  ON_NODES,
  ON_CELLS,
};

static constexpr SUPPORT_TYPE str2SupportType(std::string_view support)
{
  using enum SUPPORT_TYPE;
  if (support == "on_nodes")
    return ON_NODES;
  else if (support == "on_cells")
    return ON_CELLS;
  fmt::println(stderr, "support type {} not recognized", std::string{support});
  std::abort();

  return NONE;
}

enum struct NATURE_TYPE : uint8_t
{
  NONE = 0u,
  INTENSIVE_MAXIMUM,
  EXTENSIVE_MAXIMUM,
  INTENSIVE_CONSERVATION,
  EXTENSIVE_CONSERVATION,
};

static constexpr NATURE_TYPE str2NatureType(std::string_view nature)
{
  using enum NATURE_TYPE;
  if (nature == "intensive_maximum")
    return INTENSIVE_MAXIMUM;
  else if (nature == "extensive_maximum")
    return EXTENSIVE_MAXIMUM;
  else if (nature == "intensive_conservation")
    return INTENSIVE_CONSERVATION;
  else if (nature == "extensive_conservation")
    return EXTENSIVE_CONSERVATION;
  fmt::println(stderr, "nature type {} not recognized", std::string{nature});
  std::abort();

  return NONE;
}

struct FieldCoupling
{
  FieldCoupling() = default;
  FieldCoupling(COUPLING_TYPE type): type_{type} {}
  virtual ~FieldCoupling() = default;

  virtual size_t size() const noexcept = 0;
  virtual double const * dataPtr() const = 0;
  virtual double operator[](size_t k) const = 0;
  virtual double at(size_t k) const = 0;
  // virtual std::vector<double> getData() = 0;
  virtual void init(
      std::string_view name,
      MeshCoupling const * mesh,
      SUPPORT_TYPE support,
      NATURE_TYPE nature) = 0;
  virtual void initIO(std::filesystem::path const & prefix)
  {
    prefix_ = prefix;
    std::filesystem::create_directories(prefix_);
  }
  virtual void setValues(std::span<const double> data, uint const dim = 1u) = 0;
  virtual void setValues(double value, uint size, uint const dim = 1u) = 0;
  virtual void printVTK(double time, uint iter) = 0;

  static std::unique_ptr<FieldCoupling> build(COUPLING_TYPE type);
  static std::unique_ptr<FieldCoupling> build(std::string_view type);

  COUPLING_TYPE type_ = COUPLING_TYPE::NONE;
  std::string name_ = "";
  MeshCoupling const * mesh_;
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
      MeshCoupling const * mesh,
      SUPPORT_TYPE support,
      NATURE_TYPE nature) override
  {
    name_ = name;
    mesh_ = mesh;
    support_ = support;
    nature_ = nature;
  }
  virtual void initIO(std::filesystem::path const & /*prefix*/) override {}
  virtual void setValues(std::span<const double> data, uint const /*dim*/ = 1u) override
  {
    data_.resize(data.size());
    std::copy(data.begin(), data.end(), data_.begin());
  }
  virtual void setValues(double value, uint size, uint const /*dim*/ = 1u) override
  {
    data_.resize(size, value);
  }

  virtual void printVTK(double /*time*/, uint /*iter*/) override
  {
    if (messageVTK_)
    {
      fmt::println(stderr, "VTK print is not implemented in FieldSimple.");
      messageVTK_ = false;
    }
  }

  std::vector<double> data_;
  SUPPORT_TYPE support_;
  NATURE_TYPE nature_;
  bool messageVTK_ = true;
};

} // namespace cocoa
