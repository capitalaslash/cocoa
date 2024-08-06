#pragma once

#include <filesystem>
#include <string_view>

#include "med_mesh.hpp"

namespace MEDCoupling
{
struct MEDCouplingFieldDouble;
}

struct MEDField
{
  struct Frame
  {
    uint iter;
    double time;
  };

  MEDField() = default;
  ~MEDField();

  void init(std::string_view name, MEDMesh & mesh);
  void initIO(std::string_view filename);
  void setValues(std::vector<double> const & data);

  void printVTK(double time, uint iter);
  void printPVD();

  MEDCoupling::MEDCouplingFieldDouble * fieldPtr_;
  bool inited_ = false;
  std::filesystem::path filename_ = "tmp";
  std::vector<Frame> frames;
};
