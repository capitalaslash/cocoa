#include "med_field.hpp"

#include <string_view>

#include <fmt/ostream.h>

#include <MEDCouplingFieldDouble.hxx>

MEDField::~MEDField()
{
  if (inited_)
  {
    fieldPtr_->decrRef();
  }
}

void MEDField::init(std::string_view name, MEDMesh & mesh)
{
  fieldPtr_ = MEDCoupling::MEDCouplingFieldDouble::New(
      MEDCoupling::ON_NODES, MEDCoupling::ONE_TIME);
  inited_ = true;

  fieldPtr_->setName(std::string{name});
  fieldPtr_->setTimeUnit("s");
  fieldPtr_->setTime(0.0, 0, -1);
  fieldPtr_->setMesh(mesh.meshPtr_);
}

void MEDField::initIO(std::string_view filename) { filename_ = filename; }

void MEDField::setValues(std::vector<double> const & data)
{
  MEDCoupling::DataArrayDouble * array = MEDCoupling::DataArrayDouble::New();
  array->alloc(data.size(), 1);
  // array->fillWithValue(8.);
  std::copy(data.data(), data.data() + data.size(), array->getPointer());
  fieldPtr_->setArray(array);
  array->decrRef();
}

void MEDField::printVTK(double time, uint iter)
{
  auto const filename = filename_.string() + std::to_string(iter);
  frames.push_back(Frame{iter, time});
  fieldPtr_->setTime(time, iter, -1);
  fieldPtr_->writeVTK(filename, false);
  printPVD();
}

void MEDField::printPVD()
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