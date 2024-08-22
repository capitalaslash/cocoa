#include "coupling_simple.hpp"

// fmt
#include <fmt/ranges.h>

// local
#include "problem_fd1d.hpp"

void CouplingSimple::setup(Problem * pSrc, Problem * pTgt)
{
  // CouplingSimple works only with ProblemFD1D for now
  auto derSrc = dynamic_cast<ProblemFD1D *>(pSrc);
  assert(derSrc != nullptr);
  auto derTgt = dynamic_cast<ProblemFD1D *>(pTgt);
  assert(derTgt != nullptr);
  pSrc_ = pSrc;
  pTgt_ = pTgt;

  m_.init(derTgt->n_ - 1, derSrc->n_ - 1);
  for (uint row = 0; row < derTgt->n_ - 1; row++)
  {
    auto const startTgt = derTgt->start_ + row * derTgt->h_;
    auto const endTgt = startTgt + derTgt->h_;
    for (uint clm = 0; clm < derSrc->n_ - 1; clm++)
    {
      auto const startSrc = derSrc->start_ + clm * derSrc->h_;
      auto const endSrc = startSrc + derSrc->h_;
      if (((startSrc >= startTgt - 1.e-12) && startSrc < (endTgt + 1.e-12)) ||
          ((endSrc >= startTgt - 1.e-12) && (endSrc < endTgt + 1.e-12)))
      {
        double const start = std::max(startTgt, startSrc);
        double const end = std::min(endTgt, endSrc);
        double const frac = (end - start) / derSrc->h_;
        m_(row, clm) = frac;
      }
    }
  }
  fmt::print("m: {} x {}\n{:::.2e}\n", m_.data.size(), m_.data[0].size(), m_.data);

  // std::vector<double> fieldSrc(derSrc->n_ - 1, 1.0);
  // std::vector<double> fieldTgt(derTgt->n_ - 1, 0.0);
  // for (uint i = 0; i < fieldTgt.size(); i++)
  // {
  //   for (uint j = 0; j < fieldSrc.size(); j++)
  //   {
  //     fieldTgt[i] += m_(i, j) * fieldSrc[j];
  //   }
  // }
  // fmt::print("fieldTgt: {::.1e}\n", fieldTgt);
}

void CouplingSimple::project(
    std::string_view fieldNameSrc, std::string_view fieldNameTgt)
{
  auto derSrc = dynamic_cast<ProblemFD1D *>(pSrc_);
  auto derTgt = dynamic_cast<ProblemFD1D *>(pTgt_);

  // convert P1 to P0
  FieldCoupling * p1Src = derSrc->getField(fieldNameSrc);
  std::vector<double> p0Src(p1Src->size() - 1);
  for (uint k = 0; k < p0Src.size(); k++)
  {
    p0Src[k] = 0.5 * (*(p1Src->dataPtr() + k) + *(p1Src->dataPtr() + k + 1));
  }

  // interpolate
  FieldCoupling * p1Tgt = derTgt->getField(fieldNameTgt);
  std::vector<double> p0Tgt(p1Tgt->size() - 1);
  for (uint i = 0; i < p0Tgt.size(); i++)
  {
    for (uint j = 0; j < p0Src.size(); j++)
    {
      p0Tgt[i] += m_(i, j) * p0Src[j];
    }
  }

  // convert back P0 to P1
  std::vector<double> p1TgtData(p1Tgt->size());
  p1TgtData[0] = p0Tgt[0];
  p1TgtData[p1Tgt->size() - 1] = p0Tgt[p0Tgt.size() - 1];
  for (uint k = 1; k < p1Tgt->size() - 1; k++)
  {
    p1TgtData[k] = 0.5 * (p0Tgt[k - 1] + p0Tgt[k]);
  }
  p1Tgt->setValues(p1TgtData);
}
