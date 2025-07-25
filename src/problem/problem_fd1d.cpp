#include "problem/problem_fd1d.hpp"

// std
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <variant>

// libfmt
#include <fmt/core.h>
#include <fmt/ranges.h>

// local
#include "coupling/mesh_coupling.hpp"
#include "enums.hpp"
#include "la.hpp"
#include "problem/fdutils.hpp"

namespace cocoa
{

ProblemFD1D::ProblemFD1D(): Problem{PROBLEM_TYPE::FD1D}
{
  // register default assemblies
  assemblies_.emplace(EQN_TYPE::HEAT, [](ProblemFD1D * p) { p->assemblyHeat(); });
  assemblies_.emplace(
      EQN_TYPE::HEAT_COUPLED, [](ProblemFD1D * p) { p->assemblyHeatCoupled(); });
  assemblies_.emplace(EQN_TYPE::HEAT_OC, [](ProblemFD1D * p) { p->assemblyHeatOC(); });
}

ProblemFD1D::~ProblemFD1D()
{
  // erase possibly added assembly
  assemblies_.erase(EQN_TYPE::CUSTOM);
}

void ProblemFD1D::setup(Problem::ConfigList_T const & configs)
{
  // default values
  name_ = "empty";
  // mesh
  double start = 0.0;
  double end = 1.0;
  uint nElems = 10u;
  // fields
  nVars_ = 1u;
  varNames_ = {"u"};
  std::vector<double> uInit(nVars_, 0.0);
  std::vector<double> qValue(nVars_, 1.0);
  // bcs
  bcs_.resize(nVars_);
  // time
  time = 0.0;
  finalTime_ = 1.0;
  dt_ = 0.1;
  // linear algebra
  maxIters_ = 1000u;
  tol_ = 1.e-6;
  uint nnz = 3u;

  // read configuration from file
  // TODO: string conversion is required for python bindings, manage other variant types
  // as errors
  auto const & configFileVariant = configs.at("config_file");
  auto const configFile =
      std::holds_alternative<std::filesystem::path>(configFileVariant)
          ? std::get<std::filesystem::path>(configFileVariant)
          : std::filesystem::path{std::get<std::string>(configFileVariant)};
  std::ifstream in(configFile, std::ios::in);
  if (!in)
  {
    fmt::println(stderr, "configuration file {} not found!", configFile.string());
    std::abort();
  }
  std::string buffer;
  while (std::getline(in, buffer, '\n'))
  {
    std::istringstream bufferStream{buffer};
    std::string token;
    while (std::getline(bufferStream, token, ' '))
    {
      if (token[0] == '#')
      {
        // this is a comment, consume whole line
        while (bufferStream)
          bufferStream >> token;
      }
      else if (token == "name:")
        bufferStream >> name_;
      else if (token == "debug:")
        bufferStream >> debug_;
      // mesh
      else if (token == "start:")
        bufferStream >> start;
      else if (token == "end:")
        bufferStream >> end;
      else if (token == "n_elems:")
        bufferStream >> nElems;
      // fields
      else if (token == "n_vars:")
      {
        bufferStream >> nVars_;
        varNames_.resize(nVars_);
        uInit.resize(nVars_);
        qValue.resize(nVars_);
        bcs_.resize(nVars_);
      }
      else if (token == "var_names:")
        for (uint v = 0u; v < nVars_; v++)
          bufferStream >> varNames_[v];
      else if (token == "initial_value:")
        for (uint v = 0u; v < nVars_; v++)
          bufferStream >> uInit[v];
      else if (token == "q:")
      {
        fields_.emplace("q", VectorFD());
        for (uint v = 0u; v < nVars_; v++)
          bufferStream >> qValue[v];
      }
      else if (token == "additional_fields:")
      {
        while (bufferStream)
        {
          std::string name = "";
          bufferStream >> name;
          if (name == "")
            break;
          fields_.emplace(name, VectorFD());
        }
        fmt::print("additional fields: ");
        for (auto const & [k, v]: fields_)
          fmt::print("{} ", k);
        fmt::println("");
      }
      // params
      else if (token == "params:")
      {
        while (bufferStream)
        {
          std::string name = "";
          bufferStream >> name;
          if (name == "")
            break;
          bufferStream >> token;
          FD_PARAM_TYPE type = str2FDParamType(token);
          switch (type)
          {
          case FD_PARAM_TYPE::INTEGER:
          {
            uint value;
            bufferStream >> value;
            params_.set(name, value);
            break;
          }
          case FD_PARAM_TYPE::SCALAR:
          {
            double value;
            bufferStream >> value;
            params_.set(name, value);
            break;
          }
          case FD_PARAM_TYPE::VECTOR:
          {
            uint size;
            bufferStream >> size;
            std::vector<double> value(size);
            for (uint k = 0u; k < size; k++)
              bufferStream >> value[k];
            params_.set(name, value);
            break;
          }
          default:
            fmt::println(stderr, "param type for {} not recognized", name);
            std::abort();
          }
        }
      }
      // time
      else if (token == "start_time:")
        bufferStream >> time;
      else if (token == "final_time:")
        bufferStream >> finalTime_;
      else if (token == "dt:")
        bufferStream >> dt_;
      // assembly
      else if (token == "assembly_name:")
      {
        bufferStream >> token;
        eqnType_ = str2eqn(token);
      }
      // la
      else if (token == "solver_type:")
      {
        bufferStream >> token;
        solverType_ = str2fdsolver(token);
      }
      else if (token == "max_iters:")
        bufferStream >> maxIters_;
      else if (token == "tol:")
        bufferStream >> tol_;
      else if (token == "nnz:")
        bufferStream >> nnz;
      // bcs
      else if (token == "bc_left:")
      {
        for (uint v = 0u; v < nVars_; v++)
        {
          bufferStream >> token;
          double value;
          bufferStream >> value;
          bcs_[v].left().init(FD_BC_SIDE::LEFT, str2FDBCType(token), value);
        }
      }
      else if (token == "bc_right:")
      {
        for (uint v = 0u; v < nVars_; v++)
        {
          bufferStream >> token;
          double value;
          bufferStream >> value;
          bcs_[v].right().init(FD_BC_SIDE::RIGHT, str2FDBCType(token), value);
        }
      }
      // io
      else if (token == "print_step:")
        bufferStream >> printStep_;
      else if (token == "output_prefix:")
        bufferStream >> outputPrefix_;
      else if (token == "clean_output:")
        bufferStream >> cleanOutput_;
      else
      {
        fmt::println(stderr, "key {} invalid", token);
        bufferStream >> token;
      }
    }
  }
  fmt::println("{} - equation type: {}", name_, eqn2str(eqnType_));
  assert(eqnType_ == EQN_TYPE::NONE || assemblies_.contains(eqnType_));

  fmt::println("parameters: {}", params_);

  // mesh
  mesh_.init({start}, {end}, {nElems});

  // fields
  u_.resize(mesh_.nPts() * nVars_);
  uOld_.resize(mesh_.nPts() * nVars_);
  for (uint v = 0u; v < nVars_; v++)
  {
    u_.setRange(0u + v * mesh_.nPts(), mesh_.nPts() + v * mesh_.nPts(), uInit[v]);
    uOld_.setRange(0u + v * mesh_.nPts(), mesh_.nPts() + v * mesh_.nPts(), uInit[v]);
  }
  for (auto & [name, data]: fields_)
  {
    data.resize(mesh_.nPts() * nVars_);
    if (name == "q")
      for (uint v = 0u; v < nVars_; v++)
        fields_.at("q").setRange(
            0u + v * mesh_.nPts(), mesh_.nPts() + v * mesh_.nPts(), qValue[v]);
  }

  // linear algebra
  m_.init(mesh_.nPts() * nVars_, nnz);
  rhs_.resize(mesh_.nPts() * nVars_);

  // io
  initOutput();
}

std::unique_ptr<MeshCoupling> ProblemFD1D::initMeshCoupling(
    COUPLING_TYPE type, COUPLING_SCOPE scope, Marker marker, std::string_view bdName)
{
  if (scope == COUPLING_SCOPE::VOLUME)
  {
    // coords format: x_0, y_0, z_0, x_1, ...
    std::vector<double> coords(mesh_.nPts() * 3);
    for (uint k = 0; k < mesh_.nPts(); k++)
    {
      coords[3 * k] = mesh_.pt({k})[0];
      coords[3 * k + 1] = 0.0;
      coords[3 * k + 2] = 0.0;
    }

    // conn format: elem0_numpts, id_0, id_1, ..., elem1_numpts, ...
    auto const nElems = mesh_.nElems();
    std::vector<uint> conn(nElems * 3);
    for (uint k = 0; k < nElems; k++)
    {
      conn[3 * k] = MEDCellTypeToIKCell(MED_CELL_TYPE::LINE2);
      conn[3 * k + 1] = k;
      conn[3 * k + 2] = k + 1;
    }

    // offsets format: sum_0^k elemk_numpts + 1,
    std::vector<uint> offsets(nElems + 1);
    offsets[0] = 0;
    for (uint k = 0; k < nElems; k++)
    {
      offsets[k + 1] = offsets[k] + 3;
    }

    auto meshCoupling = MeshCoupling::build(type);
    meshCoupling->init(
        "mesh_fd1d",
        COUPLING_SCOPE::VOLUME,
        markerNotSet,
        "",
        1u,
        1u,
        coords,
        conn,
        offsets);
    return meshCoupling;
  }
  else if (scope == COUPLING_SCOPE::BOUNDARY)
  {
    fmt::println(stderr, "boundary coupling not yet implemented!");
    std::abort();
  }
  else
  {
    std::abort();
  }
}

std::unique_ptr<FieldCoupling> ProblemFD1D::initFieldCoupling(
    COUPLING_TYPE type, std::string_view name, MeshCoupling const * mesh)
{
  auto field = FieldCoupling::build(type);
  // check if the coupling field is a variable
  for (auto v = 0u; v < nVars_; v++)
    if (varNames_[v] == name)
    {
      field->init(name, mesh, SUPPORT_TYPE::ON_NODES, NATURE_TYPE::INTENSIVE_MAXIMUM);
      auto start = u_.data() + v * mesh_.nPts();
      field->setValues({start, start + mesh_.nPts()}, 1u);
      field->initIO(outputPrefix_);
      dataPtr_.emplace(name, start);
      return field;
    }
  // check if the coupling field is an additional field
  if (fields_.contains(std::string{name}))
  {
    field->init(name, mesh, SUPPORT_TYPE::ON_NODES, NATURE_TYPE::INTENSIVE_MAXIMUM);
    field->setValues(fields_.at(std::string{name}).data_, 1u);
    field->initIO(outputPrefix_);
    dataPtr_.emplace(name, fields_.at(std::string{name}).data());
    return field;
  }

  fmt::println(stderr, "coupling field {} not found", name);
  std::abort();
  return field;
}

void ProblemFD1D::setFieldData(FieldCoupling * field)
{
  // check if the field maps a variable
  for (auto v = 0u; v < nVars_; v++)
  {
    if (varNames_[v] == field->name_)
    {
      auto const start = u_.data() + v * mesh_.nPts();
      field->setValues({start, start + mesh_.nPts()});
      return;
    }
  }
  // otherwise, it must be stored in the additional fields
  field->setValues(fields_.at(field->name_).data_);
}

void ProblemFD1D::getFieldData(FieldCoupling const & field)
{
  // check if the field maps a variable
  for (auto v = 0u; v < nVars_; v++)
  {
    if (varNames_[v] == field.name_)
    {
      auto const start = u_.data() + v * mesh_.nPts();
      std::copy(field.dataPtr(), field.dataPtr() + field.size(), start);
      return;
    }
  }
  // otherwise, it must be stored in the additional fields
  std::copy(
      field.dataPtr(), field.dataPtr() + field.size(), fields_.at(field.name_).data());
}

void ProblemFD1D::initOutput()
{
  if (cleanOutput_ && std::filesystem::exists(outputPrefix_))
    for (const auto & entry: std::filesystem::directory_iterator(outputPrefix_))
      std::filesystem::remove_all(entry.path());
  std::filesystem::create_directories(outputPrefix_);
}

bool ProblemFD1D::run() const { return time < finalTime_; }

void ProblemFD1D::advance()
{
  if (time + dt_ < finalTime_ - 1.e-6)
  {
    time += dt_;
  }
  else
  {
    dt_ = finalTime_ - time;
    time = finalTime_;
  }
  it++;
}

uint ProblemFD1D::solve()
{
  fmt::println("\n===");
  fmt::println("{}, time = {:.6e}, dt = {:.6e}", name_, time, dt_);

  // update
  uOld_ = u_;

  // assembly
  assemblies_.at(eqnType_)(this);

  for (uint v = 0u; v < nVars_; v++)
  {
    // bc start
    uint const idLeft = 0u + v * mesh_.nPts();
    switch (bcs_[v].left().type)
    {
    case FD_BC_TYPE::DIRICHLET:
    {
      m_.clearRow(idLeft);
      m_.add(idLeft, idLeft, 1.0);
      rhs_.set(idLeft, bcs_[v].left().values[0]);
      break;
    }
    case FD_BC_TYPE::NEUMANN:
    {
      // sign: incoming flux is positive
      // (u_1 - u_-1) / 2h = A
      // u_-1 = u_1 - 2 h A
      // u_1 part implemented in assembly
      rhs_.add(
          idLeft,
          -2.0 * mesh_.h_[0] * bcs_[v].left().values[0] *
              bcs_[v].left().ghostValues[0]);
      break;
    }
    default:
    {
      fmt::println(stderr, "no bc left specified!");
      std::abort();
    }
    }

    // bc end
    uint const idEnd = mesh_.nPts() - 1 + v * mesh_.nPts();
    switch (bcs_[v].right().type)
    {
    case FD_BC_TYPE::DIRICHLET:
    {
      m_.clearRow(idEnd);
      m_.add(idEnd, idEnd, 1.0);
      rhs_.set(idEnd, bcs_[v].right().values[0]);
      break;
    }
    case FD_BC_TYPE::NEUMANN:
    {
      // sign: incoming flux is positive
      // (u_n-2 - u_n) / 2h = A
      // u_n = u_n-2 - 2 h A
      // u_n-2 part implemented in assembly
      rhs_.add(
          idEnd,
          -2.0 * mesh_.h_[0] * bcs_[v].right().values[0] *
              bcs_[v].right().ghostValues[0]);
      break;
    }
    default:
    {
      fmt::println(stderr, "no bc end specified!");
      std::abort();
    }
    }
  }

  m_.close();
  if (debug_)
    m_.print_sparsity_pattern("fd1d_mat.dat");

  // solve
  auto const [numIters, residual] =
      solvers_.at(solverType_)(m_, rhs_, u_, tol_, maxIters_);
  fmt::print("num iters: {:4d}, ", numIters);
  double const rhsNorm = std::sqrt(rhs_.norm2sq() * mesh_.h_[0]);
  fmt::println("relative residual: {:.8e}", residual / rhsNorm);
  if (std::isnan(residual))
  {
    fmt::println("non-scaled residual: {:.8e}", residual);
    std::abort();
  }

  if (debug_)
  {
    fmt::println("matrix: {}", m_);
    fmt::println("rhs: {}", rhs_);
    fmt::println("sol: {}", u_);
  }

  // clean up
  m_.clear();
  rhs_.zero();

  return numIters;
}

void ProblemFD1D::assemblyHeat()
{
  // eqn:
  // du/dt - alpha d^2u/dx^2 = q
  // discretization:
  // u_m / dt - alpha (u_l - 2 * u_m + u_r) / h^2 =
  // uold_m / dt + q_m
  // grouping and multiplying by dt:
  // (1 + 2 alpha dt / h^2) * u_m
  // - (alpha dt / h^2) * u_l
  // - (alpha dt / h^2) * u_r
  // = uold_m + q_m dt

  auto const alpha = params_.get<FD_PARAM_TYPE::VECTOR>("alpha");
  assert(alpha.size() == nVars_);
  auto const & q = fields_.at("q");

  for (uint v = 0u; v < nVars_; v++)
    for (uint k = 0u; k < mesh_.nPts(); k++)
    {
      uint const id = k + v * mesh_.nPts();
      auto const iH2 = 1.0 / (mesh_.h_[0] * mesh_.h_[0]);

      // diagonal
      m_.add(id, id, 1.0 + 2.0 * alpha[v] * dt_ * iH2);

      // left
      auto const valueLeft = -alpha[v] * dt_ * iH2;
      if (k > 0u)
        m_.add(id, id - 1, valueLeft);
      else
      {
        m_.add(id, id + 1, valueLeft);
        bcs_[v].left().ghostValues.set(0, valueLeft);
      }

      // right
      auto const valueRight = -alpha[v] * dt_ * iH2;
      if (k < mesh_.nPts() - 1)
        m_.add(id, id + 1, valueRight);
      else
      {
        m_.add(id, id - 1, valueRight);
        bcs_[v].right().ghostValues.set(0, valueRight);
      }

      // rhs
      rhs_.set(id, uOld_[id] + q[id] * dt_);
    }
  m_.close();
}

void ProblemFD1D::assemblyHeatCoupled()
{
  assert(nVars_ == 1u);
  auto const alpha = params_.get<FD_PARAM_TYPE::SCALAR>("alpha");
  auto const & q = fields_.at("q");

  // TODO: use VectorFD for uExt
  // std::vector<double> uExt(n_, 2.0);
  VectorFD const & uExt = fields_["Tcfd"];

  double const kAmpli = 10.;
  auto const iH2 = 1.0 / (mesh_.h_[0] * mesh_.h_[0]);
  for (uint k = 0u; k < mesh_.nPts(); k++)
  {
    // diagonal
    m_.add(
        k,
        k,
        1. / dt_               // time
            + alpha * 2. * iH2 // diffusion
            + kAmpli           // feedback control
    );

    // left
    auto const valueLeft = -alpha * iH2; // diffusion
    if (k > 0u)
      m_.add(k, k - 1, valueLeft);
    else
    {
      m_.add(k, k + 1, valueLeft);
      bcs_[0].left().ghostValues.set(0, valueLeft);
    }

    // right
    auto const valueRight = -alpha * iH2; // diffusion
    if (k < mesh_.nPts() - 1)
      m_.add(k, k + 1, valueRight);
    else
    {
      m_.add(k, k - 1, valueRight);
      bcs_[0].right().ghostValues.set(0, valueRight);
    }

    // rhs
    rhs_.set(
        k,
        uOld_[k] / dt_         // time
            + q[k]             // source
            + kAmpli * uExt[k] // feedback control
    );
  }
  m_.close();
}

void ProblemFD1D::assemblyHeatOC()
{
  auto const alpha = params_.get<FD_PARAM_TYPE::SCALAR>("alpha");
  auto const beta = params_.get<FD_PARAM_TYPE::SCALAR>("beta");
  auto const tempTarget = params_.get<FD_PARAM_TYPE::SCALAR>("tempTarget");
  auto const targetLeft = params_.get<FD_PARAM_TYPE::SCALAR>("targetLeft");
  auto const targetRight = params_.get<FD_PARAM_TYPE::SCALAR>("targetRight");

  for (uint k = 0u; k < mesh_.nPts(); k++)
  {
    auto const x = mesh_.pt({k})[0];
    auto const h = mesh_.h_[0];
    uint const idF = k;
    uint const idFLeft = (k != 0) ? idF - 1 : idF + 1;
    uint const idFRight = (k != mesh_.nPts() - 1) ? idF + 1 : idF - 1;
    uint const idA = k + mesh_.nPts();
    uint const idALeft = (k != 0) ? idA - 1 : idA + 1;
    uint const idARight = (k != mesh_.nPts() - 1) ? idA + 1 : idA - 1;

    // forward problem
    m_.add(idF, idF, 1.0 + 2.0 * alpha * dt_ / (h * h));
    m_.add(idF, idFLeft, -alpha * dt_ / (h * h));
    m_.add(idF, idFRight, -alpha * dt_ / (h * h));
    m_.add(idF, idA, dt_ / beta);
    rhs_.set(idF, uOld_[idF]);

    // adjoint problem
    m_.add(idA, idA, 1.0 + 2.0 * alpha * dt_ / (h * h));
    m_.add(idA, idALeft, -alpha * dt_ / (h * h));
    m_.add(idA, idARight, -alpha * dt_ / (h * h));
    if (x > targetLeft && x < targetRight)
    {
      m_.add(idA, idF, -dt_);
      rhs_.set(idA, -tempTarget * dt_);
    }
  }
  m_.close();
}

void ProblemFD1D::print()
{
  if (it % printStep_ == 0)
  {
    for (uint v = 0u; v < nVars_; v++)
    {
      auto const filename =
          fmt::format("{}.{}.csv", (outputPrefix_ / varNames_[v]).string(), it);
      std::FILE * out = std::fopen(filename.c_str(), "w");
      fmt::println(out, "# problemfd1d csv, mesh size: {}", mesh_.n_[0] - 1);
      fmt::print(out, "{}, ", varNames_[v]);
      for (auto const & [fieldName, _]: fields_)
        fmt::print(out, "{}, ", fieldName);
      fmt::println(out, "x, y, z");
      for (uint k = 0u; k < mesh_.nPts(); k++)
      {
        uint const id = k + v * mesh_.nPts();
        fmt::print(out, "{:.6e},", u_[id]);
        for (auto const & [_, field]: fields_)
          fmt::print(out, "{:.6e}, ", field[k]);
        fmt::println(out, "{:.6e}, 0.0, 0.0", mesh_.pt({k})[0]);
      }
      std::fclose(out);
    }
  }
}

Marker ProblemFD1D::findRegion(std::string_view name)
{
  if (name == "left")
    return 0u;
  else if (name == "right")
    return 1u;
  else
  {
    fmt::println("region {} not recognized", name);
    std::abort();
  }
  return markerNotSet;
}

std::unordered_map<FD_SOLVER_TYPE, Solver_T<ProblemFD1D::Matrix_T>>
    ProblemFD1D::solvers_ = {
        // {FD_SOLVER_TYPE::NONE,
        //  [](MatrixCSR const & m,
        //     VectorFD const & b,
        //     VectorFD & x,
        //     double const residual,
        //     uint const maxIters) {
        //    return SolverInfo{0u, 0.0};
        //  }},
        {FD_SOLVER_TYPE::JACOBI, &solveJacobi<ProblemFD1D::Matrix_T>},
        {FD_SOLVER_TYPE::GAUSS_SEIDEL, &solveGaussSeidel<ProblemFD1D::Matrix_T>},
        {FD_SOLVER_TYPE::CG, &solveConjugateGradient<ProblemFD1D::Matrix_T>},
        {FD_SOLVER_TYPE::BICGSTAB, &solveBiCGStab<ProblemFD1D::Matrix_T>},
        // {FD_SOLVER_TYPE::TRIDIAG,
        //  [](MatrixTriDiag const & m,
        //     VectorFD const & b,
        //     VectorFD & x,
        //     double const,
        //     uint const) { return solveTriDiag(m, b, x); }},
        // {FD_SOLVER_TYPE::VANKA1D, &solveVanka1D},
};

} // namespace cocoa
