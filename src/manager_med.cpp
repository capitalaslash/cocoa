#include "med_manager.hpp"

#include <fmt/ranges.h>

#include "med_mesh.hpp"
#include "problem.hpp"

void MEDManager::setup(Problem * src, Problem * tgt)
{
  problemSrc_ = src;
  problemTgt_ = tgt;
  remapper.setPrecision(1.e-12);
  remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
  remapper.prepare(
      problemSrc_->meshMED_.meshPtr_, problemTgt_->meshMED_.meshPtr_, "P1P1");
  auto const matrix = remapper.getCrudeMatrix();
  fmt::print("remapper matrix:\n{}\n", matrix);
}

void MEDManager::project(std::string_view srcName, std::string_view tgtName)
{
  MEDCoupling::MEDCouplingFieldDouble * tgtField = remapper.transferField(
      problemSrc_->fieldsMED_.at(std::string{srcName}).fieldPtr_, 0.0);
  problemTgt_->fieldsMED_.at(std::string{tgtName}).fieldPtr_ = tgtField;
}
