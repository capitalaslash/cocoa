#include "plugins.hpp"

// pybind11
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

// local
#include "coupling/coupling_med.hpp"
#include "coupling/coupling_simple.hpp"
#include "enums.hpp"
#include "problem/problem.hpp"
#include "problem/problem_fd1d.hpp"
#include "problem/problem_fd2d.hpp"
#include "problem/problem_oforg.hpp"
#include "problem/problem_proxpde.hpp"

namespace py = pybind11;

using namespace py::literals;

PYBIND11_MODULE(pycocoa, m)
{
  // enums =============================================================
  py::enum_<COUPLING_TYPE>(m, "COUPLING_TYPE")
      .value("none", COUPLING_TYPE::NONE)
      .value("medcoupling", COUPLING_TYPE::MEDCOUPLING)
      .value("ofm2m", COUPLING_TYPE::OFM2M)
      .value("simple", COUPLING_TYPE::SIMPLE);

  // problems ==========================================================
  py::class_<Problem>(m, "Problem")
      .def_readwrite("coupling_type", &Problem::couplingType_)
      .def(
          "setup",
          [](Problem * p, py::kwargs const & kwargs)
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

  py::class_<ProblemFD1D, Problem>(m, "ProblemFD1D").def(py::init<>());

  py::class_<ProblemFD2D, Problem>(m, "ProblemFD2D").def(py::init<>());

#ifdef COCOA_ENABLE_OFORG
  py::class_<ProblemOForg, Problem>(m, "ProblemOForg").def(py::init<>());
#endif

#ifdef COCOA_ENABLE_PROXPDE
  py::class_<ProblemProXPDEHeat, Problem>(m, "ProblemProXPDEHeat")
      .def(py::init<>())
      .def(
          "set_source",
          [](ProblemProXPDEHeat * p,
             std::function<double(proxpde::Vec3 const &)> const & source)
          { p->q_ << source; });

  py::class_<ProblemProXPDENS, Problem>(m, "ProblemProXPDENS").def(py::init<>());
#endif

  // coupling ==========================================================
  py::class_<CouplingManager>(m, "CouplingManager")
      .def("setup", &CouplingManager::setup, "problem_src"_a, "problem_tgt"_a)
      .def(
          "project",
          py::overload_cast<std::string_view, std::string_view>(
              &CouplingManager::project),
          "name_src"_a,
          "name_tgt"_a)
      .def(
          "project",
          py::overload_cast<std::string_view>(&CouplingManager::project),
          "name"_a);

  py::class_<CouplingSimple, CouplingManager>(m, "CouplingSimple").def(py::init<>());

#ifdef COCOA_ENABLE_MEDCOUPLING
  py::class_<CouplingMED, CouplingManager>(m, "CouplingMED").def(py::init<>());
#endif
}
