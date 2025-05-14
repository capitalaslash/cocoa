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
#include "coupling/field_coupling.hpp"
#include "coupling/field_med.hpp"
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

using namespace cocoa;

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
      .value("heat_coupled", EQN_TYPE::HEAT_COUPLED)
      .value("heat_oc", EQN_TYPE::HEAT_OC)
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
      .def("get_mesh_coupling", &Problem::getMesh)
      .def_readwrite("coupling_type", &Problem::couplingType_)
      .def_readwrite("time", &Problem::time)
      .def_readwrite("it", &Problem::it)
      .def_readwrite("print_step", &Problem::printStep_)
      .def_readwrite("debug", &Problem::debug_);

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
            p->initOutput();
          })
      .def(
          "set_custom_assembly",
          [](ProblemFD1D * p, ProblemFD1D::Assembly_T const & f)
          {
            auto const [_, success] = p->assemblies_.emplace(EQN_TYPE::CUSTOM, f);
            p->eqnType_ = EQN_TYPE::CUSTOM;
            return success;
          })
      .def_static(
          "set_custom_solver",
          [](Solver_T<ProblemFD1D::Matrix_T> const & f) -> bool
          {
            auto const [ptr, success] =
                ProblemFD1D::solvers_.emplace(FD_SOLVER_TYPE::CUSTOM, f);
            return success;
          })
      .def_static(
          "release_custom_solver",
          []() -> bool { return ProblemFD1D::solvers_.erase(FD_SOLVER_TYPE::CUSTOM); })
      .def_readwrite("name", &ProblemFD1D::name_)
      .def_readwrite("mesh", &ProblemFD1D::mesh_)
      .def_readwrite("n_vars", &ProblemFD1D::nVars_)
      .def_readwrite("var_names", &ProblemFD1D::varNames_)
      .def_readwrite("u", &ProblemFD1D::u_)
      .def_readwrite("u_old", &ProblemFD1D::uOld_)
      .def_readwrite("fields", &ProblemFD1D::fields_)
      .def_readwrite("params", &ProblemFD1D::params_)
      .def_readwrite("final_time", &ProblemFD1D::finalTime_)
      .def_readwrite("dt", &ProblemFD1D::dt_)
      .def_readwrite("m", &ProblemFD1D::m_)
      .def_readwrite("rhs", &ProblemFD1D::rhs_)
      .def_readwrite("solver_type", &ProblemFD1D::solverType_)
      .def_readwrite("max_iters", &ProblemFD1D::maxIters_)
      .def_readwrite("tol", &ProblemFD1D::tol_)
      .def_readwrite("eqn_type", &ProblemFD1D::eqnType_)
      .def_readwrite("bcs", &ProblemFD1D::bcs_)
      .def_readwrite("clean_output", &ProblemFD1D::cleanOutput_)
      .def_readwrite("output_prefix", &ProblemFD1D::outputPrefix_)
      .def_readwrite("name_ext", &ProblemFD1D::nameExt_)
      .def_readwrite("assemblies", &ProblemFD1D::assemblies_);

  // ProblemFD2D =======================================================
  py::class_<ProblemFD2D, Problem>(m, "ProblemFD2D")
      .def(py::init<>())
      .def("init_mesh_coupling", &ProblemFD2D::initMeshCoupling)
      .def("init_field_coupling", &ProblemFD2D::initFieldCoupling)
      .def(
          "create_output_dir",
          [](ProblemFD2D * p, std::string_view path)
          { std::filesystem::create_directories(path); })
      .def(
          "setup_io",
          [](ProblemFD2D * p, std::string_view path)
          {
            p->outputPrefix_ = path;
            p->initMeshCoupling();
            p->initFieldCoupling();
            p->initOutput();
          })
      .def(
          "set_custom_assembly",
          [](ProblemFD2D * p, ProblemFD2D::Assembly_T const & f)
          {
            auto const [_, success] = p->assemblies_.emplace(EQN_TYPE::CUSTOM, f);
            p->eqnType_ = EQN_TYPE::CUSTOM;
            return success;
          })
      .def_readwrite("name", &ProblemFD2D::name_)
      .def_readwrite("mesh", &ProblemFD2D::mesh_)
      .def_readwrite("n_vars", &ProblemFD2D::nVars_)
      .def_readwrite("var_names", &ProblemFD2D::varNames_)
      .def_readwrite("u", &ProblemFD2D::u_)
      .def_readwrite("u_old", &ProblemFD2D::uOld_)
      .def_readwrite("fields", &ProblemFD2D::fields_)
      .def_readwrite("params", &ProblemFD2D::params_)
      .def_readwrite("q", &ProblemFD2D::q_)
      .def_readwrite("final_time", &ProblemFD2D::finalTime_)
      .def_readwrite("dt", &ProblemFD2D::dt_)
      .def_readwrite("m", &ProblemFD2D::m_)
      .def_readwrite("rhs", &ProblemFD2D::rhs_)
      .def_readwrite("solver_type", &ProblemFD2D::solverType_)
      .def_readwrite("max_iters", &ProblemFD2D::maxIters_)
      .def_readwrite("tol", &ProblemFD2D::tol_)
      .def_readwrite("eqn_type", &ProblemFD2D::eqnType_)
      .def_readwrite("bcs", &ProblemFD2D::bcs_)
      .def_readwrite("output_prefix", &ProblemFD2D::outputPrefix_)
      .def_readwrite("name_ext", &ProblemFD2D::nameExt_)
      .def_readwrite("assemblies", &ProblemFD2D::assemblies_);

  // ProblemOForg ======================================================
#ifdef COCOA_ENABLE_OFORG
  py::class_<ProblemOForg, Problem>(m, "ProblemOForg").def(py::init<>());
#endif

  // ProblemProXPDE ====================================================
#ifdef COCOA_ENABLE_PROXPDE
  py::class_<ProblemProXPDEHeat, Problem>(m, "ProblemProXPDEHeat")
      .def(py::init<>())
      .def(
          "set_source",
          [](ProblemProXPDEHeat * p,
             std::function<double(proxpde::Vec3 const &)> const & source)
          {
            p->fieldsP0_.emplace("q", proxpde::FEVar{"q", p->feSpaceP0_});
            p->fieldsP0_.at("q") << source;
          });

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

  py::class_<FieldCoupling>(m, "FieldCoupling")
      // .def("init", &FieldCoupling::init, "name"_a, "mesh"_a, "support"_a)
      .def(
          "init",
          [](FieldCoupling * f,
             std::string_view name,
             Problem * problem,
             std::string_view support)
          {
            if (support == "on_nodes")
              f->init(name, problem->getMesh(), SUPPORT_TYPE::ON_NODES);
            else if (support == "on_cells")
              f->init(name, problem->getMesh(), SUPPORT_TYPE::ON_CELLS);
            else
            {
              fmt::print(stderr, "support type {} not recognized\n", support);
              std::abort();
            }
          })
      .def("at", &FieldCoupling::at, "pos"_a)
      .def(
          "set_values",
          [](FieldCoupling * f, std::vector<double> & data, uint const dim)
          { f->setValues(data, dim); },
          "data"_a,
          "dim"_a = 1u)
      .def(
          "set_values",
          py::overload_cast<double, uint, uint const>(&FieldCoupling::setValues),
          "value"_a,
          "size"_a,
          "dim"_a = 1u)
      .def("init_io", &FieldCoupling::initIO, "prefix"_a)
      .def("print_vtk", &FieldCoupling::printVTK, "time"_a, "iter"_a);

#ifdef COCOA_ENABLE_MEDCOUPLING
  py::class_<FieldMED, FieldCoupling>(m, "FieldMED").def(py::init<>());
#endif

  // FD utils ==========================================================
  py::enum_<FD_SOLVER_TYPE>(m, "FD_SOLVER_TYPE")
      .value("jacobi", FD_SOLVER_TYPE::JACOBI)
      .value("gauss_seidel", FD_SOLVER_TYPE::GAUSS_SEIDEL)
      .value("cg", FD_SOLVER_TYPE::CG)
      .value("custom", FD_SOLVER_TYPE::CUSTOM);

  py::enum_<FD_BC_SIDE>(m, "FD_BC_SIDE")
      .value("left", FD_BC_SIDE::LEFT)
      .value("right", FD_BC_SIDE::RIGHT)
      .value("bottom", FD_BC_SIDE::BOTTOM)
      .value("top", FD_BC_SIDE::TOP)
      .value("front", FD_BC_SIDE::FRONT)
      .value("back", FD_BC_SIDE::BACK);

  py::enum_<FD_BC_TYPE>(m, "FD_BC_TYPE")
      .value("none", FD_BC_TYPE::NONE)
      .value("dirichlet", FD_BC_TYPE::DIRICHLET)
      .value("neumann", FD_BC_TYPE::NEUMANN);

  py::class_<FDBC>(m, "FDBC")
      .def(py::init<FD_BC_SIDE, FD_BC_TYPE, VectorFD>(), "side"_a, "type"_a, "values"_a)
      .def(
          py::init<FD_BC_SIDE, FD_BC_TYPE, double, size_t>(),
          "side"_a,
          "type"_a,
          "value"_a,
          "size"_a = 1u)
      .def_readwrite("values", &FDBC::values)
      .def_readwrite("ghost_values", &FDBC::ghostValues);

  py::class_<FDBCList1D>(m, "FDBCList1D")
      .def(py::init<>())
      .def(py::init<FDBC const &, FDBC const &>(), "left"_a, "right"_a)
      .def_property("left", &FDBCList1D::left, &FDBCList1D::left)
      .def_property("right", &FDBCList1D::right, &FDBCList1D::right)
      .def_readwrite("data", &FDBCList1D::data_);

  py::class_<FDBCList2D>(m, "FDBCList2D")
      .def(py::init<>())
      .def(
          py::init<FDBC const &, FDBC const &, FDBC const &, FDBC const &>(),
          "left"_a,
          "right"_a,
          "bottom"_a,
          "top"_a)
      .def_property("bottom", &FDBCList2D::bottom, &FDBCList2D::bottom)
      .def_property("right", &FDBCList2D::right, &FDBCList2D::right)
      .def_property("top", &FDBCList2D::top, &FDBCList2D::top)
      .def_property("left", &FDBCList2D::left, &FDBCList2D::left)
      .def_readwrite("data", &FDBCList2D::data_);

  py::class_<MeshFD1D>(m, "MeshFD1D")
      .def(py::init<>())
      .def(
          py::init<
              MeshFD1D::Real_T const &,
              MeshFD1D::Real_T const &,
              MeshFD1D::Int_T const &>(),
          "start"_a,
          "end"_a,
          "n"_a)
      .def("init", &MeshFD1D::init, "start"_a, "end"_a, "n"_a)
      // .def("pt", py::overload_cast<MeshFD1D::Int_T const &>(&MeshFD1D::pt),
      // "ijk"_a)
      .def("pt", &MeshFD1D::pt, "ijk"_a)
      .def_property_readonly("n_pts", &MeshFD1D::nPts)
      .def_readwrite("h", &MeshFD1D::h_);

  py::class_<MeshFD2D>(m, "MeshFD2D")
      .def(py::init<>())
      .def(
          py::init<
              MeshFD2D::Real_T const &,
              MeshFD2D::Real_T const &,
              MeshFD2D::Int_T const &>(),
          "start"_a,
          "end"_a,
          "n"_a)
      .def("init", &MeshFD2D::init, "start"_a, "end"_a, "n"_a)
      // .def("pt", py::overload_cast<MeshFD2D::Int_T const &>(&MeshFD2D::pt),
      // "ijk"_a)
      .def("pt", &MeshFD2D::pt, "ijk"_a)

      .def_property_readonly("n_pts", &MeshFD2D::nPts)
      .def_readwrite("h", &MeshFD2D::h_)
      .def_readwrite("n", &MeshFD2D::n_);

  // ParamsFD ==========================================================
  py::enum_<FD_PARAM_TYPE>(m, "FD_PARAM_TYPE")
      .value("integer", FD_PARAM_TYPE::INTEGER)
      .value("scalar", FD_PARAM_TYPE::SCALAR)
      .value("vector", FD_PARAM_TYPE::VECTOR);

  py::class_<ParamsFD>(m, "ParamsFD")
      .def(py::init<>())
      .def("set_integer", &ParamsFD::set<uint>, "name"_a, "value"_a)
      .def("set_scalar", &ParamsFD::set<double>, "name"_a, "value"_a)
      .def("set_vector", &ParamsFD::set<std::vector<double>>, "name"_a, "value"_a)
      .def("get_integer", &ParamsFD::get<FD_PARAM_TYPE::INTEGER>, "name"_a)
      .def("get_scalar", &ParamsFD::get<FD_PARAM_TYPE::SCALAR>, "name"_a)
      .def("get_vector", &ParamsFD::get<FD_PARAM_TYPE::VECTOR>, "name"_a);

  // linear algebra ====================================================
  py::class_<VectorFD>(m, "VectorFD", py::buffer_protocol())
      .def(py::init<size_t, double>(), "size"_a, "value"_a = 0.0)
      .def(py::init<std::vector<double>>())
      .def("resize", &VectorFD::resize, "size"_a, "value"_a = 0.0)
      .def("__getitem__", &VectorFD::operator[], "index"_a)
      .def("add", &VectorFD::add, "index"_a, "value"_a)
      .def("set", &VectorFD::set, "index"_a, "value"_a)
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
          })
      .def_readwrite("data", &VectorFD::data_);

  py::class_<MatrixTriDiag>(m, "MatrixTriDiag")
      .def(py::init<size_t>())
      .def("init", &MatrixTriDiag::init, "n"_a)
      .def("add", &MatrixTriDiag::add, "row"_a, "clm"_a, "value"_a)
      .def("close", &MatrixTriDiag::close);

  py::class_<MatrixCSR::Entry>(m, "CSREntry")
      // .def("__iter__", []() {})
      .def_readwrite("clm", &MatrixCSR::Entry::clm)
      .def_readwrite("value", &MatrixCSR::Entry::value);

  py::class_<MatrixCSR>(m, "MatrixCSR")
      .def(py::init<size_t, size_t>())
      .def("init", &MatrixCSR::init, "n"_a, "nnz"_a)
      .def("add", &MatrixCSR::add, "row"_a, "clm"_a, "value"_a)
      .def("close", &MatrixCSR::close)
      .def_readwrite("data", &MatrixCSR::data_);

  py::class_<SolverInfo>(m, "SolverInfo")
      .def(py::init<uint, double>(), "n_iters"_a, "residual"_a)
      .def_readwrite("n_iters", &SolverInfo::nIters)
      .def_readwrite("residual", &SolverInfo::residual);
}
