#include "med_mesh.hpp"

#include <MEDCouplingUMesh.hxx>

MEDMesh::~MEDMesh()
{
  if (inited_)
  {
    meshPtr_->decrRef();
  }
}

void MEDMesh::init(
    std::string_view name,
    uint const dim,
    std::vector<double> const & coords,
    std::vector<mcIdType> conn,
    std::vector<mcIdType> offsets)
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
}

void MEDMesh::printVTK() { meshPtr_->writeVTK(meshPtr_->getName(), false); }
