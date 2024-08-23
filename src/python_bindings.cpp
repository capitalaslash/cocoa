// pybind11
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// local
#include "coupling_med.hpp"
#include "coupling_simple.hpp"
#include "enums.hpp"
#include "problem.hpp"
#include "problem_fd1d.hpp"
#include "problem_proxpde.hpp"

using namespace pybind11::literals;

PYBIND11_MODULE(pycocoa, m)
{
  pybind11::class_<std::filesystem::path>(m, "Path").def(pybind11::init<std::string>());
  pybind11::implicitly_convertible<std::string, std::filesystem::path>();

  pybind11::enum_<COUPLING_TYPE>(m, "COUPLING_TYPE")
      .value("none", COUPLING_TYPE::NONE)
      .value("medcoupling", COUPLING_TYPE::MEDCOUPLING)
      .value("ofm2m", COUPLING_TYPE::OFM2M)
      .value("simple", COUPLING_TYPE::SIMPLE);

  pybind11::class_<Problem>(m, "Problem")
      .def_readwrite("coupling_type", &Problem::couplingType_)
      .def(
          "setup",
          [](Problem * p, pybind11::kwargs const & kwargs)
          {
            using ParamList_T = Problem::ParamList_T;
            using Key_T = ParamList_T::key_type;
            // using Value_T = ParamList_T::value_type;
            ParamList_T configFiles;
            for (auto kw: kwargs)
            {
              configFiles[kw.first.cast<Key_T>()] = kw.second.cast<std::string>();
            }
            p->setup(configFiles);
          })
      .def("run", &Problem::run)
      .def("advance", &Problem::advance)
      .def("solve", &Problem::solve)
      .def("print", &Problem::print)
      .def("getField", &Problem::getField, "name"_a)
      .def("setField", &Problem::setField, "name"_a, "field"_a);

  m.def("setFD1DAssemblies", &setFD1DAssemblies, "Initializa FD1D assenbly routines.");

  pybind11::class_<ProblemFD1D, Problem>(m, "ProblemFD1D").def(pybind11::init<>());

  pybind11::class_<ProblemProXPDE, Problem>(m, "ProblemProXPDE")
      .def(pybind11::init<>());

  pybind11::class_<CouplingManager>(m, "CouplingManager")
      .def("setup", &CouplingManager::setup, "problem_src"_a, "problem_tgt"_a)
      .def("project", &CouplingManager::project, "name_src"_a, "name_tgt"_a);

  pybind11::class_<CouplingSimple, CouplingManager>(m, "CouplingSimple")
      .def(pybind11::init<>());

  pybind11::class_<CouplingMED, CouplingManager>(m, "CouplingMED")
      .def(pybind11::init<>());
}
