#include <string>
#include <string_view>
#include <vector>

#include "med_field.hpp"
#include "problem.hpp"

struct ProblemFD1D: public Problem
{
  struct Matrix
  {
    Matrix(uint n): diag(n), diagUp(n), diagDown(n) {}
    std::vector<double> diag;
    std::vector<double> diagUp;
    std::vector<double> diagDown;
  };

  ProblemFD1D() = default;
  ~ProblemFD1D() = default;

  void setup(Problem::ParamList_T const & params) override;
  bool run() override;
  void advance() override;
  void solve() override;
  void print() override;

  void initMeshMED(std::filesystem::path const & fileName);
  void initFieldMED(std::filesystem::path const & fileName);

  MEDField getField(std::string_view name) override { return MEDField{}; }
  void setField(std::string_view name, MEDField const & field) override {}

  std::vector<double> points_;
  MEDMesh meshMED_;
  std::vector<double> u_;
  MEDField uMED_;
  std::vector<double> uOld_;
  std::vector<double> q_;
  double diff_;
  double finalTime_;
  double dt_;
  double bcStart_;
  double bcEnd_;
  std::string outFile_ = "fd1d";
};
