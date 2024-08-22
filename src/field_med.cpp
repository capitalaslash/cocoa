#include "field_med.hpp"

// std
#include <string_view>

// fmt
#include <fmt/ostream.h>

// medcoupling
#include <MEDCouplingFieldDouble.hxx>

// local
#include "mesh_med.hpp"

FieldMED::~FieldMED()
{
  if (inited_)
  {
    fieldPtr_->decrRef();
  }
}

// std::vector<double> FieldMED::getData() override { return data_; }

void FieldMED::init(std::string_view name, MeshCoupling * mesh)
{
  fieldPtr_ = MEDCoupling::MEDCouplingFieldDouble::New(
      MEDCoupling::ON_NODES, MEDCoupling::ONE_TIME);
  inited_ = true;

  // TODO: allow more field natures
  fieldPtr_->setNature(MEDCoupling::IntensiveMaximum);
  fieldPtr_->setName(std::string{name});
  fieldPtr_->setTimeUnit("s");
  fieldPtr_->setTime(0.0, 0, -1);
  fieldPtr_->setMesh(dynamic_cast<MeshMED *>(mesh)->meshPtr_);
}

void FieldMED::initIO(std::string_view filename) { filename_ = filename; }

void FieldMED::setValues(std::vector<double> const & data)
{
  MEDCoupling::DataArrayDouble * array = MEDCoupling::DataArrayDouble::New();
  array->alloc(data.size(), 1);
  std::copy(data.data(), data.data() + data.size(), array->getPointer());
  fieldPtr_->setArray(array);
  array->decrRef();
}

void FieldMED::setValues(double value, uint size)
{
  MEDCoupling::DataArrayDouble * array = MEDCoupling::DataArrayDouble::New();
  array->alloc(size, 1);
  array->fillWithValue(value);
  fieldPtr_->setArray(array);
  array->decrRef();
}

void FieldMED::printVTK(double time, uint iter)
{
  auto const filename = filename_.string() + std::to_string(iter);
  frames.push_back(Frame{iter, time});
  fieldPtr_->setTime(time, iter, -1);
  fieldPtr_->writeVTK(filename, false);
  printPVD();
}

void FieldMED::printPVD() const
{
  std::FILE * out = std::fopen((filename_.string() + "pvd").c_str(), "w");
  fmt::print(
      out,
      "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" "
      "header_type=\"UInt64\">\n");
  fmt::print(out, "  <Collection>\n");
  for (auto const & frame: frames)
  {
    fmt::print(
        out,
        "    <DataSet timestep=\"{:.6e}\" part=\"0\" file=\"{}{}.vtu\"/>\n",
        frame.time,
        filename_.filename().string(),
        frame.iter);
  }
  fmt::print(out, "  </Collection>\n");
  fmt::print(out, "</VTKFile>\n");
  std::fclose(out);
}