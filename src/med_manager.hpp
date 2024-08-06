#pragma once

#include <vector>

#include "med_field.hpp"
#include "med_mesh.hpp"

struct Problem;

struct MEDManager
{
  Problem * problemSrc_;
  Problem * problemTgt_;
};
