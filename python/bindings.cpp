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
#include "la.hpp"
#include "problem/fdutils.hpp"
#include "problem/problem.hpp"
#include "problem/problem_fd1d.hpp"
#include "problem/problem_fd2d.hpp"
#include "problem/problem_oforg.hpp"
#include "problem/problem_proxpde.hpp"

namespace py = pybind11;

using namespace py::literals;

PYBIND11_MODULE(cocoa, m)
{
  // enums =============================================================
  py::enum_<COUPLING_TYPE>(m, "COUPLING_TYPE")
      .value("none", COUPLING_TYPE::NONE)
      .value("medcoupling", COUPLING_TYPE::MEDCOUPLING)
      .value("ofm2m", COUPLING_TYPE::OFM2M)
      .value("simple", COUPLING_TYPE::SIMPLE);

  py::enum_<EQN_TYPE>(m, "EQN_TYPE")
      .value("none", EQN_TYPE::NONE)
      .value("heat", EQN_TYPE::HEAT)
      .value("custom", EQN_TYPE::CUSTOM);

  // Problem ===========================================================
  py::class_<Problem>(m, "Problem")
      .def(
          "setup",
          [](Problem * p, py::kwargs const & kwargs)
          {
            using ConfigList_T = Problem::ConfigList_T;
            using Key_T = ConfigList_T::key_type;
            // using Value_T = ConfigList_T::value_type;
            ConfigList_T configFiles;
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
      .def("get_field", &Problem::getField, "name"_a)
      .def("set_field", &Problem::setField, "name"_a, "field"_a)
      .def_readwrite("coupling_type", &Problem::couplingType_)
      .def_readwrite("time", &Problem::time);

  // ProblemFD1D =======================================================
  py::class_<ProblemFD1D, Problem>(m, "ProblemFD1D")
      .def(py::init<>())
      .def("init_mesh_coupling", &ProblemFD1D::initMeshCoupling)
      .def("init_field_coupling", &ProblemFD1D::initFieldCoupling)
      .def(
          "create_output_dir",
          [](ProblemFD1D * p, std::string_view path)
          { std::filesystem::create_directories(path); })
      .def(
          "setup_io",
          [](ProblemFD1D * p, std::string_view path)
          {
            p->outputPrefix_ = path;
            p->initMeshCoupling();
            p->initFieldCoupling();
            std::filesystem::create_directories(path);
          })
      .def(
          "set_custom_assembly",
          [](ProblemFD1D * p, ProblemFD1D::Assembly_T const & f)
          {
            auto const [_, success] = p->assemblies_.emplace(EQN_TYPE::CUSTOM, f);
            p->eqnType_ = EQN_TYPE::CUSTOM;
            return success;
          })
      .def_readwrite("name", &ProblemFD1D::name_)
      .def_readwrite("start", &ProblemFD1D::start_)
      .def_readwrite("h", &ProblemFD1D::h_)
      .def_readwrite("n", &ProblemFD1D::n_)
      .def_readwrite("var_name", &ProblemFD1D::varName_)
      .def_readwrite("u", &ProblemFD1D::u_)
      .def_readwrite("u_old", &ProblemFD1D::uOld_)
      .def_readwrite("q", &ProblemFD1D::q_)
      .def_readwrite("alpha", &ProblemFD1D::alpha_)
      .def_readwrite("final_time", &ProblemFD1D::finalTime_)
      .def_readwrite("dt", &ProblemFD1D::dt_)
      .def_readwrite("m", &ProblemFD1D::m_)
      .def_readwrite("rhs", &ProblemFD1D::rhs_)
      .def_readwrite("solver_type", &ProblemFD1D::solverType_)
      .def_readwrite("eqn_type", &ProblemFD1D::eqnType_)
      .def_readwrite("bc_start", &ProblemFD1D::bcStart_)
      .def_readwrite("bc_end", &ProblemFD1D::bcEnd_)
      .def_readwrite("output_prefix", &ProblemFD1D::outputPrefix_)
      .def_readwrite("name_ext", &ProblemFD1D::nameExt_)
      .def_readwrite("assemblies", &ProblemFD1D::assemblies_);

  // other derived problems ============================================
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

  // FD utils ==========================================================
  py::enum_<FD_BC_TYPE>(m, "FD_BC_TYPE")
      .value("none", FD_BC_TYPE::NONE)
      .value("dirichlet", FD_BC_TYPE::DIRICHLET)
      .value("neumann", FD_BC_TYPE::NEUMANN);

  py::class_<FDBC>(m, "FDBC")
      .def(py::init<FD_BC_TYPE, std::vector<double>>())
      .def_readwrite("ghost_values", &FDBC::ghostValues);

  // linear algebra ====================================================
  py::class_<VectorFD>(m, "VectorFD", py::buffer_protocol())
      .def(py::init<size_t>())
      .def_buffer(
          [](VectorFD & v)
          {
            return py::buffer_info(
                v.data(),
                sizeof(double),
                py::format_descriptor<double>::format(),
                1U,
                {v.size()},
                {sizeof(double)});
          });

  py::class_<MatrixTriDiag>(m, "MatrixTriDiag")
      .def(py::init<size_t>())
      .def("init", &MatrixTriDiag::init, "n"_a)
      .def("add", &MatrixTriDiag::add, "row"_a, "clm"_a, "value"_a)
      .def("close", &MatrixTriDiag::close);

  py::class_<MatrixCSR>(m, "MatrixCSR")
      .def(py::init<size_t>())
      .def("init", &MatrixCSR::init, "n"_a)
      .def("add", &MatrixCSR::add, "row"_a, "clm"_a, "value"_a)
      .def("close", &MatrixCSR::close);
}
