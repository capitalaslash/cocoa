#include "coupling_med.hpp"

// std
#include <cassert>

// fmt
#include <fmt/ranges.h>

// local
#include "coupling/field_med.hpp"
#include "coupling/mesh_med.hpp"
#include "problem/problem.hpp"

namespace cocoa
{

void CouplingMED::setup(Problem * src, Problem * tgt)
{
  problemSrc_ = src;
  problemTgt_ = tgt;

  remapper.setPrecision(1.e-12);
  remapper.setIntersectionType(interpType_);
  auto meshSrc = dynamic_cast<MeshMED *>(problemSrc_->meshCoupling_.get());
  auto meshTgt = dynamic_cast<MeshMED *>(problemTgt_->meshCoupling_.get());
  remapper.prepare(meshSrc->meshPtr_, meshTgt->meshPtr_, "P1P1");

  auto const matrix = remapper.getCrudeMatrix();
  // fmt::println("remapper matrix:\n{}", matrix);

  auto mask = std::vector<double>(matrix.size());
  for (uint k = 0; k < matrix.size(); k++)
  {
    auto sum = 0.0;
    for (auto & [k, v]: matrix[k])
    {
      sum += v;
    }
    mask[k] = (sum > 1.e-12) ? 1.0 : 0.0;
  }

  auto [kvPair, success] = problemTgt_->fieldsCoupling_.emplace("mask", new FieldMED);
  assert(success);
  kvPair->second->init("mask", meshTgt, SUPPORT_TYPE::ON_NODES);
  kvPair->second->setValues(mask);
}

void CouplingMED::project(std::string_view srcName, std::string_view tgtName)
{
  auto srcField = dynamic_cast<FieldMED *>(problemSrc_->getField(srcName));
  auto tgtField = dynamic_cast<FieldMED *>(problemTgt_->getField(tgtName));
  // MEDCoupling::MEDCouplingFieldDouble * tmp =
  //     remapper.transferField(srcField->fieldPtr_, 0.0);
  // tgtField->fieldPtr_ = tmp;
  tgtField->fieldPtr_ = remapper.transferField(srcField->fieldPtr_, 0.0);
}

} // namespace cocoa
