#pragma once

// std
#include <filesystem>
#include <string>
#include <vector>

// fmt
#include <fmt/core.h>

// local
#include "enums.hpp"
#include "mesh_coupling.hpp"

struct FieldCoupling
{
  FieldCoupling() = default;
  explicit FieldCoupling(COUPLING_TYPE type): type_{type} {}
  virtual ~FieldCoupling() = default;

  virtual size_t size() const noexcept = 0;
  virtual double const * dataPtr() = 0;
  // virtual std::vector<double> getData() = 0;
  virtual void init(std::string_view name, MeshCoupling * mesh) = 0;
  virtual void initIO(std::string_view filename) = 0;
  virtual void setValues(std::vector<double> const & data) = 0;
  virtual void setValues(double value, uint size) = 0;
  virtual void printVTK(double time, uint iter) = 0;

  COUPLING_TYPE type_ = COUPLING_TYPE::NONE;
  std::string name_ = "";
  std::filesystem::path filename_ = "./tmp";
};

struct FieldSimple: public FieldCoupling
{
  FieldSimple(): FieldCoupling{COUPLING_TYPE::SIMPLE} {}
  ~FieldSimple() = default;

  size_t size() const noexcept override { return data_.size(); }
  double const * dataPtr() override { return data_.data(); }
  // std::vector<double> getData() override { return data_; }
  virtual void init(std::string_view name, MeshCoupling * mesh) override {}
  virtual void initIO(std::string_view filename) override {}
  virtual void setValues(std::vector<double> const & data) override { data_ = data; }
  virtual void setValues(double value, uint size) override
  {
    data_.resize(size, value);
  }
  virtual void printVTK(double time, uint iter) override
  {
    if (messageVTK_)
    {
      fmt::print(stderr, "VTK print is not implemented in CouplingSimple.\n");
      messageVTK_ = false;
    }
  }

  std::vector<double> data_;
  bool messageVTK_ = true;
};
