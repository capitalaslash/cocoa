#pragma once

// std
#include <filesystem>
#include <string>
#include <vector>

// local
#include "mesh_coupling.hpp"

struct FieldCoupling
{
  FieldCoupling() = default;
  virtual ~FieldCoupling() = default;

  virtual size_t size() const noexcept { return 0; }
  virtual double * dataPtr() { return nullptr; }
  virtual void init(std::string_view name, MeshCoupling * mesh) {}
  virtual void initIO(std::string_view filename) {}
  virtual void setValues(std::vector<double> const & data) {}
  virtual void setValues(double value, uint size) {}

  std::string name_;
  std::filesystem::path filename_ = "tmp";
};

struct FieldSimple: public FieldCoupling
{
  FieldSimple() = default;
  ~FieldSimple() = default;

  size_t size() const noexcept override { return data_.size(); }
  double * dataPtr() override { return data_.data(); }
  virtual void init(std::string_view name, MeshCoupling * mesh) override {}
  virtual void initIO(std::string_view filename) override {}
  virtual void setValues(std::vector<double> const & data) override { data_ = data; }
  virtual void setValues(double value, uint size) override
  {
    data_.resize(size, value);
  }

  std::vector<double> data_;
};