#include "problem/problem_oforg.hpp"

// std
#include <array>
#include <cassert>
#include <fstream>
#include <unordered_map>

// openfoam.org
#include <IOobject.H>
#include <volFieldsFwd.H>
// #include <ConstantField.H>  // only in openfoam.com

// local
#include "enums.hpp"

namespace cocoa
{

void ProblemOForg::setup(Problem::ConfigList_T const & configs)
{
  // read case dir
  // TODO: string conversion is required for python bindings, manage other variant types
  // as errors
  auto const & caseDirVariant = configs.at("case_dir");
  prefix_ = std::holds_alternative<std::filesystem::path>(caseDirVariant)
                ? std::get<std::filesystem::path>(caseDirVariant)
                : std::filesystem::path{std::get<std::string>(caseDirVariant)};
  int argc = 3;
  char ** argv = new char *[3];
  argv[0] = (char *)"app_oforg";
  argv[1] = (char *)"-case";
  argv[2] = new char[prefix_.string().size()];
  std::strcpy(argv[2], prefix_.string().data());

  Foam::argList::addOption("solver", "name", "Solver name");

  Foam::argList args(argc, argv);

  if (!args.checkRootCase())
  {
    Foam::FatalError.exit();
  }

  Foam::Info << "Create time\n" << Foam::endl;

  runTime_.reset(new Foam::Time{Foam::Time::controlDictName, args});

  // Read the solverName from the optional solver entry in controlDict
  solverName_ = runTime_->controlDict().lookupOrDefault("solver", Foam::word::null);

  // Optionally reset the solver name from the -solver command-line argument
  args.optionReadIfPresent("solver", solverName_);

  // Check the solverName has been set
  if (solverName_ == Foam::word::null)
  {
    args.printUsage();

    FatalErrorIn(args.executable())
        << "solver not specified in the controlDict or on the command-line"
        << exit(Foam::FatalError);
  }
  else
  {
    // Load the solver library
    // TODO: add map with custom lib names
    Foam::solver::load(solverName_);
  }

  // Create the default single region mesh
  Foam::Info << "Create mesh for time = " << runTime_->name() << Foam::nl << Foam::endl;

  mesh_.reset(new Foam::fvMesh{Foam::IOobject{
      Foam::fvMesh::defaultRegion,
      runTime_->name(),
      *runTime_,
      Foam::IOobject::MUST_READ}});

  initMeshMED("mesh_oforg");

  // Instantiate the selected solver
  solverPtr_ = Foam::solver::New(solverName_, *mesh_);
  // Foam::solver & solver = solverPtr();

  // Create the outer PIMPLE loop and control structure
  pimple_.reset(new Foam::pimpleSingleRegionControl{solverPtr_->pimple});

  // Set the initial time-step
  setDeltaT(*runTime_, *solverPtr_);

  std::string outputVTK = "./output_oforg";
  if (configs.contains("config_file"))
  {
    auto const & configFileVariant = configs.at("config_file");
    auto const configFile =
        std::holds_alternative<std::filesystem::path>(configFileVariant)
            ? std::get<std::filesystem::path>(configFileVariant)
            : std::filesystem::path{std::get<std::string>(configFileVariant)};
    std::ifstream in(configFile, std::ios::in);
    if (!in)
    {
      fmt::print(stderr, "configuration file {} not found!\n", configFile.string());
      std::abort();
    }
    std::string buffer;
    while (std::getline(in, buffer, '\n'))
    {
      std::istringstream bufferStream{buffer};
      std::string token;
      while (std::getline(bufferStream, token, ' '))
      {
        if (token == "scalar_vars:")
          while (std::getline(bufferStream, token, ' '))
          {
            namesExport_.emplace_back(token, OFFIELD_TYPE::SCALAR);
          }
        else if (token == "vector_vars:")
          while (std::getline(bufferStream, token, ' '))
          {
            namesExport_.emplace_back(token, OFFIELD_TYPE::VECTOR);
          }
        else if (token == "output_vtk:")
          bufferStream >> outputVTK;
        else
        {
          fmt::print(stderr, "key {} invalid\n", token);
          bufferStream >> token;
        }
      }
    }
  }

  for (auto const & [name, type]: namesExport_)
  {
    initFieldMED(name, outputVTK + "/" + name);
    switch (type)
    {
    case OFFIELD_TYPE::SCALAR:
    {
      Foam::volScalarField const & field =
          mesh_->lookupObject<Foam::volScalarField>(name);
      setDataMED(name, field);
      break;
    }
    case OFFIELD_TYPE::VECTOR:
    {
      Foam::volVectorField const & field =
          mesh_->lookupObject<Foam::volVectorField>(name);
      setDataMED(name, field);
      break;
    }
    default:
      std::abort();
    }
    getField(name)->printVTK(0.0, 0);
  }

  // field_.reset(new Foam::volScalarField{
  //     Foam::IOobject{
  //         "field",
  //         runTime_->name(),
  //         *mesh_,
  //         Foam::IOobject::NO_READ,
  //         Foam::IOobject::NO_WRITE},
  //     *mesh_,
  //     Foam::dimensionedScalar{Foam::dimless, 2.0}});
  // forAll(mesh_->cells(), cellId) { (*field_)[cellId] = cellId; }

  Foam::Info << Foam::nl << "Starting time loop\n" << Foam::endl;
}

// TODO: add enum
static std::unordered_map<uint8_t, std::array<uint, 8U>> connOF2MED = {
    {8U, {0U, 1U, 3U, 2U, 4U, 5U, 7U, 6U}}, // HEX8
};

void ProblemOForg::initMeshMED(std::string_view meshName)
{
  // // boundary
  //   Foam::label patchId = mesh_->boundaryMesh().findPatchID("name");
  //   const Foam::polyPatch & namePolyPatch = mesh_->boundaryMesh()[patchId];
  //   Foam::labelList test(namePolyPatch.boundaryPoints());
  //   forAll(test, pos) { Foam::Info << test[pos] << Foam::endl; }

  // coords format: x_0, y_0, z_0, x_1, ...
  std::vector<double> coords(mesh_->nPoints() * 3);
  forAll(mesh_->points(), pointId)
  {
    auto const & pt = mesh_->points()[pointId];
    coords[3 * pointId + 0] = pt[0];
    coords[3 * pointId + 1] = pt[1];
    coords[3 * pointId + 2] = pt[2];
  }

  // conn format: elem0_numpts, id_0, id_1, ..., elem1_numpts, ...
  // TODO: extend for generic cells, not only hex8
  std::vector<uint> conn(mesh_->nCells() * (8 + 1));
  forAll(mesh_->cells(), cellId)
  {
    conn[(8 + 1) * cellId] = MEDCellTypeToIKCell(MED_CELL_TYPE::HEX8);
    for (uint p = 0; p < 8; p++)
    {
      conn[(8 + 1) * cellId + 1 + connOF2MED[8U][p]] = mesh_->cellPoints()[cellId][p];
    }
  }

  // offsets format: sum_0^k elemk_numpts + 1,
  std::vector<uint> offsets(mesh_->nCells() + 1);
  offsets[0] = 0;
  forAll(mesh_->cells(), cellId) { offsets[cellId + 1] = offsets[cellId] + 8 + 1; }

  couplingType_ = COUPLING_TYPE::MEDCOUPLING;
  meshCoupling_ = MeshCoupling::build(couplingType_);
  meshCoupling_->init(meshName, 3U, 3U, coords, conn, offsets);
  meshCoupling_->printVTK(prefix_ / "mesh_oforg");
}

// TODO: move to Problem since it does not require any specific knowledge of the
// specific type
void ProblemOForg::initFieldMED(std::string_view name, std::filesystem::path prefix)
{
  auto [kvPair, success] =
      fieldsCoupling_.emplace(name, FieldCoupling::build(COUPLING_TYPE::MEDCOUPLING));
  assert(success);
  kvPair->second->init(name, meshCoupling_.get(), SUPPORT_TYPE::ON_CELLS);
  kvPair->second->initIO(prefix);
}

template <typename Field>
void ProblemOForg::setDataMED(std::string_view fieldName, Field const & field)
{
  // static_assert(std::is_same_v<Field, Foam::GeometricField>);
  std::vector<double> data;
  if constexpr (std::is_same_v<typename Field::value_type, Foam::scalar>)
  {
    data.resize(field.size());
    forAll(mesh_->cells(), cellId) { data[cellId] = field[cellId]; }
    getField(fieldName)->setValues(data, 1U);
  }
  else if constexpr (std::is_same_v<typename Field::value_type, Foam::vector>)
  {
    std::vector<double> data(field.size() * 3U);
    forAll(mesh_->cells(), cellId)
    {
      data[3 * cellId + 0] = field[cellId][0];
      data[3 * cellId + 1] = field[cellId][1];
      data[3 * cellId + 2] = field[cellId][2];
    }
    getField(fieldName)->setValues(data, 3U);
  }
  else
  {
    std::abort();
  }
}

bool ProblemOForg::run() { return pimple_->run(*runTime_); }

void ProblemOForg::advance()
{
  // Update PIMPLE outer-loop parameters if changed
  pimple_->read();

  solverPtr_->preSolve();

  // Adjust the time-step according to the solver maxDeltaT
  adjustDeltaT(*runTime_, *solverPtr_);

  (*runTime_)++;

  Foam::Info << "Time = " << runTime_->userTimeName() << Foam::nl << Foam::endl;
}

uint ProblemOForg::solve()
{
  // PIMPLE corrector loop
  while (pimple_->loop())
  {
    solverPtr_->moveMesh();
    // solverPtr_->motionCorrector();
    solverPtr_->fvModels().correct();
    solverPtr_->prePredictor();
    solverPtr_->momentumPredictor();
    solverPtr_->thermophysicalPredictor();
    solverPtr_->pressureCorrector();
    solverPtr_->postCorrector();
  }

  solverPtr_->postSolve();

  time += runTime_->deltaTValue();
  it++;

  for (auto const & [name, type]: namesExport_)
  {
    switch (type)
    {
    case OFFIELD_TYPE::SCALAR:
    {
      Foam::volScalarField const & field =
          mesh_->lookupObject<Foam::volScalarField>(name);
      setDataMED(name, field);
      break;
    }
    case OFFIELD_TYPE::VECTOR:
    {
      Foam::volVectorField const & field =
          mesh_->lookupObject<Foam::volVectorField>(name);
      setDataMED(name, field);
      break;
    }
    default:
      std::abort();
    }
  }

  // TODO: extract the number of iterations
  return 0u;
}

void ProblemOForg::print()
{
  runTime_->write();

  Foam::Info << "ExecutionTime = " << runTime_->elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime_->elapsedClockTime() << " s" << Foam::nl
             << Foam::endl;

  bool writeVTK = false;
  Foam::word const writeStyle = runTime_->controlDict().lookup("writeControl");
  if (writeStyle == "timeStep")
  {
    Foam::label writeInterval =
        runTime_->controlDict().lookupOrDefault("writeInterval", 1);
    if (it % writeInterval == 0)
      writeVTK = true;
  }
  else if (writeStyle == "runTime")
  {
    double const dtWrite =
        runTime_->controlDict().lookupOrDefault("writeInterval", 1.0);
    if (lastPrint_ + dtWrite > time - 1.e-12)
    {
      writeVTK = true;
      lastPrint_ = time;
    }
  }
  else
  {
    fmt::print(stderr, "writeControl {} not supported\n", writeStyle);
  }

  if (writeVTK)
    for (auto const & [name, type]: namesExport_)
    {
      getField(name)->printVTK(time, it);
    }
}

void setDeltaT(Foam::Time & runTime, const Foam::solver & solver)
{
  if (runTime.timeIndex() == 0 &&
      runTime.controlDict().lookupOrDefault("adjustTimeStep", false) &&
      solver.transient())
  {
    // const Foam::scalar deltaT =
    //     Foam::min(solver.maxDeltaT(), runTime.functionObjects().maxDeltaT());
    const Foam::scalar deltaT = solver.maxDeltaT();

    if (deltaT < Foam::rootVGreat)
    {
      runTime.setDeltaT(Foam::min(runTime.deltaTValue(), deltaT));
    }
  }
}

void adjustDeltaT(Foam::Time & runTime, const Foam::solver & solver)
{
  // Update the time-step limited by the solver maxDeltaT
  if (runTime.controlDict().lookupOrDefault("adjustTimeStep", false) &&
      solver.transient())
  {
    // const Foam::scalar deltaT =
    //     Foam::min(solver.maxDeltaT(), runTime.functionObjects().maxDeltaT());
    const Foam::scalar deltaT = solver.maxDeltaT();

    if (deltaT < Foam::rootVGreat)
    {
      runTime.setDeltaT(
          Foam::min(Foam::solver::deltaTFactor * runTime.deltaTValue(), deltaT));
      Foam::Info << "deltaT = " << runTime.deltaTValue() << Foam::endl;
    }
  }
}

} // namespace cocoa
