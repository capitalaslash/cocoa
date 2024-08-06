#include <string>
#include <string_view>
#include <vector>

#define HAVE_MED 1

#ifdef HAVE_MED
#include "med_field.hpp"
#endif

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

  void setup(Problem::ParamList_T const & params) override;
  bool run() override;
  void advance() override;
  void solve() override;
  void print() override;

#ifdef HAVE_MED
  MEDField get_field(std::string_view name) override { return MEDField{}; }

  void set_field(std::string_view name, MEDField const & field) override {}
#endif

  std::vector<double> points_;
  std::vector<double> u_;
  std::vector<double> uOld_;
  std::vector<double> q_;
  double diff_;
  double finalTime_;
  double dt_;
  double bcStart_;
  double bcEnd_;
  std::string outFile_ = "fd1d";
};
