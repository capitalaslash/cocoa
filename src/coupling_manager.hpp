#pragma once

// std
#include <memory>
#include <string_view>

// local
#include "enums.hpp"

struct Problem;

struct CouplingManager
{
  CouplingManager() = default;
  explicit CouplingManager(COUPLING_TYPE type): type_(type) {}
  virtual ~CouplingManager() = default;

  virtual void setup(Problem * pSrc, Problem * pTgt) = 0;

  virtual void
  project(std::string_view fieldNameSrc, std::string_view fieldNameTgt) = 0;
  virtual void project(std::string_view fieldName) final
  {
    this->project(fieldName, fieldName);
  }

  static std::unique_ptr<CouplingManager> build(COUPLING_TYPE type);
  static std::unique_ptr<CouplingManager> build(std::string_view type);

  COUPLING_TYPE type_ = COUPLING_TYPE::NONE;
  Problem * pSrc_;
  Problem * pTgt_;
};
