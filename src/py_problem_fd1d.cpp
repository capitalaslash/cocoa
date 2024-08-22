#include "problem_fd1d.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(fd1d, m)
{
  pybind11::class_<std::filesystem::path>(m, "Path").def(pybind11::init<std::string>());
  pybind11::implicitly_convertible<std::string, std::filesystem::path>();

  pybind11::class_<Problem>(m, "Problem");

  pybind11::class_<ProblemFD1D, Problem>(m, "ProblemFD1D")
      .def(pybind11::init<>())
      .def("setup", &ProblemFD1D::setup)
      .def("run", &ProblemFD1D::run)
      .def("advance", &ProblemFD1D::advance)
      .def("solve", &ProblemFD1D::solve)
      .def("print", &ProblemFD1D::print)
      .def("getField", &ProblemFD1D::getField)
      .def("setField", &ProblemFD1D::setField);

  m.def("setFD1DAssemblies", &setFD1DAssemblies, "Initializa FD1D assenbly routines.");
}
