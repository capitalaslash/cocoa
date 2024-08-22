#include "coupling_med.hpp"

// fmt
#include <fmt/ranges.h>

// local
#include "field_med.hpp"
#include "mesh_med.hpp"
#include "problem.hpp"

void CouplingMED::setup(Problem * src, Problem * tgt)
{
  problemSrc_ = src;
  problemTgt_ = tgt;
  remapper.setPrecision(1.e-12);
  remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
  auto meshSrc = dynamic_cast<MeshMED *>(problemSrc_->meshCoupling_.get());
  auto meshTgt = dynamic_cast<MeshMED *>(problemTgt_->meshCoupling_.get());
  remapper.prepare(meshSrc->meshPtr_, meshTgt->meshPtr_, "P1P1");
  auto const matrix = remapper.getCrudeMatrix();
  fmt::print("remapper matrix:\n{}\n", matrix);
}

void CouplingMED::project(std::string_view srcName, std::string_view tgtName)
{
  auto srcField = dynamic_cast<FieldMED *>(
      problemSrc_->fieldsCoupling_.at(std::string{srcName}).get());
  auto tgtField = dynamic_cast<FieldMED *>(
      problemTgt_->fieldsCoupling_.at(std::string{tgtName}).get());
  MEDCoupling::MEDCouplingFieldDouble * tmp =
      remapper.transferField(srcField->fieldPtr_, 0.0);
  tgtField->fieldPtr_ = tmp;
}
