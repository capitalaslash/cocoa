#pragma once

#include <string_view>
#include <vector>

#include <MCIdType.hxx>
#include <MEDCouplingUMesh.hxx>
// missing extension!!!
#include <NormalizedGeometricTypes>

enum struct MED_CELL_TYPE : char
{
  LINE2 = INTERP_KERNEL::NORM_SEG2,
  TRIANGLE3 = INTERP_KERNEL::NORM_TRI3,
  QUAD4 = INTERP_KERNEL::NORM_QUAD4,
};

inline mcIdType MEDCellTypeToIKCell(MED_CELL_TYPE t)
{
  return static_cast<mcIdType>(t);
}

struct MeshMED: public MeshCoupling
{
  MeshMED() = default;
  ~MeshMED();

  void init(
      std::string_view name,
      uint const dim,
      std::vector<double> const & coords,
      std::vector<uint> conn,
      std::vector<uint> offsets) override;

  void printVTK();

  MEDCoupling::MEDCouplingUMesh * meshPtr_;
  bool inited_ = false;
};
