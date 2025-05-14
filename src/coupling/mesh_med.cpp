#include "coupling/mesh_med.hpp"

// medcoupling
#include <MEDCouplingUMesh.hxx>

namespace cocoa
{

MeshMED::~MeshMED()
{
  if (inited_)
  {
    meshPtr_->decrRef();
  }
}

void MeshMED::init(
    std::string_view name,
    uint const dim,
    uint const spaceDim,
    std::vector<double> const & coords,
    std::vector<uint> conn,
    std::vector<uint> offsets)
{
  meshPtr_ = MEDCoupling::MEDCouplingUMesh::New(std::string{name}, dim);
  inited_ = true;

  MEDCoupling::DataArrayIdType * nodalConn = MEDCoupling::DataArrayIdType::New();
  nodalConn->alloc(conn.size(), 1);
  std::copy(conn.data(), conn.data() + conn.size(), nodalConn->getPointer());
  MEDCoupling::DataArrayIdType * nodalConnI = MEDCoupling::DataArrayIdType::New();
  nodalConnI->alloc(offsets.size(), 1);
  std::copy(offsets.data(), offsets.data() + offsets.size(), nodalConnI->getPointer());
  meshPtr_->setConnectivity(nodalConn, nodalConnI, true);
  nodalConn->decrRef();
  nodalConnI->decrRef();

  MEDCoupling::DataArrayDouble * coordsArr = MEDCoupling::DataArrayDouble::New();
  coordsArr->alloc(coords.size() / 3, 3);
  std::copy(coords.data(), coords.data() + coords.size(), coordsArr->getPointer());
  meshPtr_->setCoords(coordsArr);
  coordsArr->decrRef();
  meshPtr_->changeSpaceDimension(spaceDim);
}

void MeshMED::printVTK(std::filesystem::path const & path)
{
  meshPtr_->writeVTK(path.string(), false);
}

} // namespace cocoa
