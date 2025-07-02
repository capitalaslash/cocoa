#include "coupling/field_med.hpp"

// std
#include <string_view>

// fmt
#include <fmt/ostream.h>

// medcoupling
#include <MEDCouplingFieldDouble.hxx>

// local
#include "coupling/mesh_med.hpp"

namespace cocoa
{

FieldMED::~FieldMED()
{
  if (inited_)
  {
    fieldPtr_->decrRef();
  }
}

// std::vector<double> FieldMED::getData() override { return data_; }

void FieldMED::init(
    std::string_view name, MeshCoupling const * mesh, SUPPORT_TYPE support)
{
  name_ = name;
  fieldPtr_ = MEDCoupling::MEDCouplingFieldDouble::New(
      supportType2MEDtype(support), MEDCoupling::ONE_TIME);
  inited_ = true;

  // TODO: allow more field natures
  fieldPtr_->setNature(MEDCoupling::IntensiveMaximum);
  fieldPtr_->setName(std::string{name});
  fieldPtr_->setTimeUnit("s");
  fieldPtr_->setTime(0.0, 0, -1);
  fieldPtr_->setMesh(dynamic_cast<MeshMED const *>(mesh)->meshPtr_);
}

void FieldMED::setValues(std::span<const double> data, uint const dim)
{
  MEDCoupling::DataArrayDouble * array = MEDCoupling::DataArrayDouble::New();
  array->alloc(data.size() / dim, dim);
  std::copy(data.data(), data.data() + data.size(), array->getPointer());
  fieldPtr_->setArray(array);
  array->decrRef();
}

void FieldMED::setValues(double value, uint size, uint const dim)
{
  MEDCoupling::DataArrayDouble * array = MEDCoupling::DataArrayDouble::New();
  array->alloc(size / dim, dim);
  array->fillWithValue(value);
  fieldPtr_->setArray(array);
  array->decrRef();
}

void FieldMED::printVTK(double time, uint iter)
{
  auto const filename = fmt::format("{}.{}", (prefix_ / name_).string(), iter);
  frames.push_back(Frame{iter, time});
  fieldPtr_->setTime(time, iter, -1);
  fieldPtr_->writeVTK(filename, /*binary*/ false);
  printPVD();
}

void FieldMED::printPVD() const
{
  auto const filenamePVD = fmt::format("{}.pvd", (prefix_ / name_).string());
  std::FILE * out = std::fopen(filenamePVD.c_str(), "w");
  fmt::println(
      out,
      "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" "
      "header_type=\"UInt64\">");
  fmt::println(out, "  <Collection>");
  for (auto const & frame: frames)
  {
    auto const filename = fmt::format("{}.{}.vtu", name_, frame.iter);
    fmt::println(
        out,
        "    <DataSet timestep=\"{:.6e}\" part=\"0\" file=\"{}\"/>",
        frame.time,
        filename);
  }
  fmt::println(out, "  </Collection>");
  fmt::println(out, "</VTKFile>");
  std::fclose(out);
}

} // namespace cocoa
