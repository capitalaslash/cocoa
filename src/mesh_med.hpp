#pragma once

// std
#include <string_view>
#include <vector>

// medcoupling
#include <MCIdType.hxx>
#include <MEDCouplingUMesh.hxx>
#include <NormalizedGeometricTypes> // missing extension!!!

// local
#include "mesh_coupling.hpp"

enum struct MED_CELL_TYPE : int8_t
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
  MeshMED(): MeshCoupling(COUPLING_TYPE::MEDCOUPLING) {}
  ~MeshMED();

  void init(
      std::string_view name,
      uint const dim,
      uint const spaceDim,
      std::vector<double> const & coords,
      std::vector<uint> conn,
      std::vector<uint> offsets) override;

  void printVTK();

  MEDCoupling::MEDCouplingUMesh * meshPtr_;
  bool inited_ = false;
};
