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

void CouplingMED::setup(CouplingInterface interfaceSrc, CouplingInterface interfaceTgt)
{
  interfaces_[0] = interfaceSrc;
  interfaces_[1] = interfaceTgt;

  remapper_.setPrecision(1.e-12);
  remapper_.setIntersectionType(interpType_);
  meshes_[0] = interfaces_[0].ptr->initMeshCoupling(COUPLING_TYPE::MEDCOUPLING);
  meshes_[1] = interfaces_[1].ptr->initMeshCoupling(COUPLING_TYPE::MEDCOUPLING);
  auto * meshSrc = dynamic_cast<MeshMED *>(meshes_[0].get());
  auto * meshTgt = dynamic_cast<MeshMED *>(meshes_[1].get());
  remapper_.prepare(meshSrc->meshPtr_, meshTgt->meshPtr_, "P1P1");

  auto const matrix = remapper_.getCrudeMatrix();
  // fmt::println(stderr, "remapper_ matrix:\n{}", matrix);

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

  mask_.init("mask", meshTgt, SUPPORT_TYPE::ON_NODES);
  mask_.setValues(mask);

  initFieldCoupling();
}

void CouplingMED::initFieldCoupling()
{
  for (auto k = 0u; k < 2u; k++)
  {
    for (auto const & fieldName: interfaces_[k].fieldNames)
    {
      auto field =
          interfaces_[k].ptr->initFieldCoupling(type_, fieldName, meshes_[k].get());
      auto [_, success] = fields_[k].emplace(fieldName, field.release());
      assert(success);
    }
  }
}

void CouplingMED::project(std::string_view nameSrc, std::string_view nameTgt)
{
  getField(nameSrc);
  auto srcField = dynamic_cast<FieldMED *>(fields_[0].at(std::string{nameSrc}).get());
  // srcField->fieldPtr_->writeVTK("./src", false);
  auto tgtField = dynamic_cast<FieldMED *>(fields_[1].at(std::string{nameTgt}).get());
  tgtField->fieldPtr_ = remapper_.transferField(srcField->fieldPtr_, /*default*/ 0.0);
  // tgtField->fieldPtr_->writeVTK("./tgt", false);
  setField(nameTgt);
}

} // namespace cocoa
