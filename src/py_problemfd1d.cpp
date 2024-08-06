#include "problem_fd1d.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(fd1d, m)
{
  py::class_<std::filesystem::path>(m, "Path").def(py::init<std::string>());
  py::implicitly_convertible<std::string, std::filesystem::path>();

  py::class_<ProblemFD1D>(m, "ProblemFD1D")
      .def(py::init<>())
      .def("setup", &ProblemFD1D::setup)
      .def("run", &ProblemFD1D::run)
      .def("advance", &ProblemFD1D::advance)
      .def("solve", &ProblemFD1D::solve)
      .def("print", &ProblemFD1D::print);
}
