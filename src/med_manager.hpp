#pragma once

#include <vector>

#include "coupling_manager.hpp"
#include "med_field.hpp"
#include "med_mesh.hpp"

struct Problem;

struct MEDManager: public CouplingManager
{
  Problem * problemSrc_;
  Problem * problemTgt_;

  void project(std::string_view fieldNameSrc, std::string_view fieldNameTgt) override;
};
