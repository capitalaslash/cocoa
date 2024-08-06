#pragma once

#include <string_view>
#include <vector>

#include <MCIdType.hxx>
#include <MEDCouplingUMesh.hxx>

struct MEDMesh
{
  MEDMesh() = default;

  ~MEDMesh();

  void init(
      std::string_view name,
      std::vector<double> const & coords,
      std::vector<mcIdType> conn,
      std::vector<mcIdType> offsets);

  void printVtk();

  MEDCoupling::MEDCouplingUMesh * meshPtr_;
  bool inited_ = false;
};
