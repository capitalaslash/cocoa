#include "coupling_simple.hpp"

// std
#include <cassert>

// fmt
#include <fmt/ranges.h>

// local
#include "enums.hpp"
#include "problem/problem_fd1d.hpp"

namespace cocoa
{

void CouplingSimple::setup(
    CouplingInterface interfaceSrc, CouplingInterface interfaceTgt)
{
  // CouplingSimple works only with ProblemFD1D for now
  auto derSrc = dynamic_cast<ProblemFD1D *>(interfaceSrc.ptr);
  assert(derSrc != nullptr);
  auto derTgt = dynamic_cast<ProblemFD1D *>(interfaceTgt.ptr);
  assert(derTgt != nullptr);
  interfaces_[0] = interfaceSrc;
  interfaces_[1] = interfaceTgt;

  meshes_[0] = interfaces_[0].ptr->initMeshCoupling(COUPLING_TYPE::SIMPLE);
  meshes_[1] = interfaces_[1].ptr->initMeshCoupling(COUPLING_TYPE::SIMPLE);
  initFieldCoupling();

  // P1P1
  m_.init(derTgt->mesh_.nPts(), derSrc->mesh_.nPts());
  for (uint row = 0u; row < derTgt->mesh_.nPts(); row++)
  {
    auto const ptTgt = derTgt->mesh_.pt({row})[0];
    for (uint clm = 0u; clm < derSrc->mesh_.nPts() - 1; clm++)
    {
      auto const startSrc = derSrc->mesh_.pt({clm})[0];
      auto const endSrc = derSrc->mesh_.pt({clm + 1})[0];
      if (ptTgt >= startSrc && ptTgt < endSrc)
      {
        m_(row, clm) = (endSrc - ptTgt) / (endSrc - startSrc);
        m_(row, clm + 1) = 1. - m_(row, clm);
        break;
      }
    }
  }

  // P0P0
  // m_.init(derTgt->n_ - 1, derSrc->n_ - 1);
  // for (uint row = 0; row < derTgt->n_ - 1; row++)
  // {
  //   auto const startTgt = derTgt->start_ + row * derTgt->h_;
  //   auto const endTgt = startTgt + derTgt->h_;
  //   for (uint clm = 0; clm < derSrc->n_ - 1; clm++)
  //   {
  //     auto const startSrc = derSrc->start_ + clm * derSrc->h_;
  //     auto const endSrc = startSrc + derSrc->h_;
  //     if (((startSrc >= startTgt - 1.e-12) && startSrc < (endTgt + 1.e-12)) ||
  //         ((endSrc >= startTgt - 1.e-12) && (endSrc < endTgt + 1.e-12)))
  //     {
  //       double const start = std::max(startTgt, startSrc);
  //       double const end = std::min(endTgt, endSrc);
  //       double const frac = (end - start) / derSrc->h_;
  //       m_(row, clm) = frac;
  //     }
  //   }
  // }

  // fmt::println("m: {} x {}\n{:::.2e}", m_.data.size(), m_.data[0].size(), m_.data);
}

void CouplingSimple::initFieldCoupling()
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

void CouplingSimple::project(
    std::string_view fieldNameSrc, std::string_view fieldNameTgt)
{
  getField(fieldNameSrc);
  // interpolate P1P1
  FieldCoupling const & p1Src = *(fields_[0].at(std::string{fieldNameSrc}));
  FieldCoupling & p1Tgt = *(fields_[1].at(std::string{fieldNameTgt}));
  std::vector<double> dataTgt(p1Tgt.size());
  for (uint i = 0; i < p1Tgt.size(); i++)
    for (uint j = 0; j < p1Src.size(); j++)
    {
      dataTgt[i] += m_(i, j) * p1Src[j];
    }
  p1Tgt.setValues(dataTgt);
  setField(fieldNameTgt);

  // // convert P1 to P0
  // FieldCoupling const & p1Src = *derSrc->getField(fieldNameSrc);
  // std::vector<double> p0Src(p1Src.size() - 1);
  // for (uint k = 0; k < p0Src.size(); k++)
  // {
  //   p0Src[k] = 0.5 * (p1Src[k] + p1Src[k + 1]);
  // }

  // // interpolate P0P0
  // FieldCoupling & p1Tgt = *derTgt->getField(fieldNameTgt);
  // std::vector<double> p0Tgt(p1Tgt->size() - 1);
  // for (uint i = 0; i < p0Tgt.size(); i++)
  // {
  //   for (uint j = 0; j < p0Src.size(); j++)
  //   {
  //     p0Tgt[i] += m_(i, j) * p0Src[j];
  //   }
  // }

  // // convert back P0 to P1
  // std::vector<double> p1TgtData(p1Tgt->size());
  // p1TgtData[0] = p0Tgt[0];
  // p1TgtData[p1Tgt->size() - 1] = p0Tgt[p0Tgt.size() - 1];
  // for (uint k = 1; k < p1Tgt.size() - 1; k++)
  // {
  //   p1TgtData[k] = 0.5 * (p0Tgt[k - 1] + p0Tgt[k]);
  // }
  // p1Tgt.setValues(p1TgtData);
}

} // namespace cocoa
