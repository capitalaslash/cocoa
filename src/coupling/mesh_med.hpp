#pragma once

#include "plugins.hpp"

#ifdef COCOA_ENABLE_MEDCOUPLING

// std
#include <string_view>
#include <vector>

// medcoupling
#include <MCIdType.hxx>
#include <MEDCouplingUMesh.hxx>
#include <NormalizedGeometricTypes> // missing extension!!!

// local
#include "coupling/mesh_coupling.hpp"

namespace cocoa
{

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

  void printVTK(std::filesystem::path const & path) override;

  MEDCoupling::MEDCouplingUMesh * meshPtr_;
  bool inited_ = false;
};

} // namespace cocoa

#endif
