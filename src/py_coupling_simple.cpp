#include "coupling_simple.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(coupling_simple, m)
{
  pybind11::class_<CouplingSimple>(m, "CouplingSimple")
      .def(pybind11::init<>())
      .def("setup", &CouplingSimple::setup)
      .def("project", &CouplingSimple::project);
}