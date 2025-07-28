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

bool ProblemOForg::argInit = true;

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
  std::strncpy(argv[2], prefix_.string().data(), prefix_.string().size());

  std::vector<std::pair<std::string, OFFIELD_TYPE>> namesExport = {};
  outputVTK_ = "./output_oforg";
  bool blockMesh = false;
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
        if (token == "scalar_vars:")
          while (std::getline(bufferStream, token, ' '))
          {
            namesExport.emplace_back(token, OFFIELD_TYPE::SCALAR);
          }
        else if (token == "vector_vars:")
          while (std::getline(bufferStream, token, ' '))
          {
            namesExport.emplace_back(token, OFFIELD_TYPE::VECTOR);
          }
        else if (token == "output_vtk:")
          bufferStream >> outputVTK_;
        else if (token == "block_mesh:")
          bufferStream >> blockMesh;
        else
        {
          fmt::println(stderr, "key {} invalid", token);
          bufferStream >> token;
        }
      }
    }
  }

  Foam::argList::addOption("solver", "name", "Solver name");

  args_.reset(new Foam::argList(
      argc, argv, /*checkArgs=*/true, /*checkOpts=*/true, /*initialize=*/argInit));
  argInit = false;

  if (blockMesh)
    runBlockMesh(*args_);

  if (!args_->checkRootCase())
  {
    Foam::FatalError.exit();
  }

  Foam::Info << "Create time\n" << Foam::endl;

  runTime_.reset(new Foam::Time{Foam::Time::controlDictName, *args_});

  // Read the solverName from the optional solver entry in controlDict
  solverName_ = runTime_->controlDict().lookupOrDefault("solver", Foam::word::null);

  // Optionally reset the solver name from the -solver command-line argument
  args_->optionReadIfPresent("solver", solverName_);

  // Check the solverName has been set
  if (solverName_ == Foam::word::null)
  {
    args_->printUsage();

    FatalErrorIn(args_->executable())
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

  // Instantiate the selected solver
  solverPtr_ = Foam::solver::New(solverName_, *mesh_);
  // Foam::solver & solver = solverPtr();

  // Create the outer PIMPLE loop and control structure
  pimple_.reset(new Foam::pimpleSingleRegionControl{solverPtr_->pimple});

  // Set the initial time-step
  setDeltaT(*runTime_, *solverPtr_);

  meshExport_ = initMeshCouplingVolume(COUPLING_TYPE::MEDCOUPLING);

  for (auto const & [name, type]: namesExport)
  {
    auto const & [it, success] = fieldsExport_.emplace(
        name, initFieldCoupling(COUPLING_TYPE::MEDCOUPLING, name, meshExport_.get()));
    assert(success);
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
static std::unordered_map<uint8_t, std::vector<uint>> connOF2MED = {
    {4u, {0u, 1u, 2u, 3u}},                 // QUAD4
    {8u, {0u, 1u, 3u, 2u, 4u, 5u, 7u, 6u}}, // HEX8
};

std::unique_ptr<MeshCoupling> ProblemOForg::initMeshCoupling(
    COUPLING_TYPE type, COUPLING_SCOPE scope, Marker marker, std::string_view bdName)
{
  // only MEDCOUPLING supported for now
  assert(type == COUPLING_TYPE::MEDCOUPLING);

  if (scope == COUPLING_SCOPE::VOLUME)
  {
    // coords format: x_0, y_0, z_0, x_1, ...
    std::vector<double> coords(mesh_->nPoints() * 3);
    for (Foam::label pointId = 0; pointId < mesh_->nPoints(); pointId++)
    {
      auto const & pt = mesh_->points()[pointId];
      coords[3 * pointId + 0] = pt[0];
      coords[3 * pointId + 1] = pt[1];
      coords[3 * pointId + 2] = pt[2];
    }

    // conn format: elem0_numpts, id_0, id_1, ..., elem1_numpts, ...
    // TODO: extend for generic cells, not only hex8
    std::vector<uint> conn(mesh_->nCells() * (8 + 1));
    for (Foam::label cellId = 0; cellId < mesh_->cells().size(); cellId++)
    {
      conn[(8 + 1) * cellId] = MEDCellTypeToIKCell(MED_CELL_TYPE::HEX8);
      for (auto p = 0u; p < 8u; p++)
        conn[(8 + 1) * cellId + 1 + connOF2MED.at(8u)[p]] =
            mesh_->cellPoints()[cellId][p];
    }

    // offsets format: sum_0^k elemk_numpts + 1,
    std::vector<uint> offsets(mesh_->nCells() + 1);
    offsets[0] = 0;
    for (Foam::label cellId = 0; cellId < mesh_->cells().size(); cellId++)
      offsets[cellId + 1] = offsets[cellId] + 8 + 1;

    auto mesh = MeshCoupling::build(type);
    mesh->init(
        "mesh_oforg",
        COUPLING_SCOPE::VOLUME,
        markerNotSet,
        "",
        3u,
        3u,
        coords,
        conn,
        offsets);
    // mesh->printVTK(outputVTK_ / "mesh_oforg");

    return mesh;
  }
  else if (scope == COUPLING_SCOPE::BOUNDARY)
  {
    // when coupling on boundary, we need a marker that identifies the patch
    assert(marker != markerNotSet);

    const Foam::polyPatch & patch = mesh_->boundaryMesh()[marker];

    // // coords format: x_0, y_0, z_0, x_1, ...
    std::vector<double> coords(patch.nPoints() * 3);
    auto const & patchPtIds = patch.meshPoints();
    std::unordered_map<Foam::label, Foam::label> globalToLocal;
    for (Foam::label localId = 0; localId < patch.nPoints(); localId++)
    {
      auto const globalId = patchPtIds[localId];
      auto const & pt = mesh_->points()[globalId];
      coords[3 * localId + 0] = pt[0];
      coords[3 * localId + 1] = pt[1];
      coords[3 * localId + 2] = pt[2];
      globalToLocal[globalId] = localId;
    }

    // conn format: elem0_numpts, id_0, id_1, ..., elem1_numpts, ...
    // TODO: extend for generic cells, not only quad4
    // using mesh->cellShapes()
    std::vector<uint> conn(patch.size() * (4 + 1));
    // auto const patchStart = patch.start();
    for (Foam::label cellId = 0; cellId < patch.size(); cellId++)
    {
      conn[(4 + 1) * cellId] = MEDCellTypeToIKCell(MED_CELL_TYPE::QUAD4);
      for (auto p = 0u; p < 4u; p++)
      {
        auto const globalId = patch[cellId][p];
        conn[(4 + 1) * cellId + 1 + connOF2MED.at(4u)[p]] = globalToLocal.at(globalId);
      }
    }

    // // offsets format: sum_0^k elemk_numpts + 1,
    std::vector<uint> offsets(patch.size() + 1);
    offsets[0] = 0;
    for (Foam::label cellId = 0; cellId < patch.size(); cellId++)
      offsets[cellId + 1] = offsets[cellId] + 4 + 1;

    auto meshBd = MeshCoupling::build(type);
    meshBd->init(
        "meshbd_oforg",
        COUPLING_SCOPE::BOUNDARY,
        marker,
        bdName,
        2u,
        3u,
        coords,
        conn,
        offsets);
    // meshBd->printVTK(outputVTK_ / "meshbd_oforg");

    return meshBd;
  }
  else
  {
    // should never get here
    std::abort();
  }
  return MeshCoupling::build(type);
}

std::unique_ptr<FieldCoupling> ProblemOForg::initFieldCoupling(
    COUPLING_TYPE type, std::string_view name, MeshCoupling const * mesh)
{
  auto field = FieldCoupling::build(type);
  field->init(name, mesh, SUPPORT_TYPE::ON_CELLS, NATURE_TYPE::INTENSIVE_MAXIMUM);
  setFieldData(field.get());
  field->initIO(outputVTK_);
  return field;
}

void ProblemOForg::setFieldData(FieldCoupling * fieldTgt)
{
  if (fieldTgt->mesh_->scope_ == COUPLING_SCOPE::VOLUME)
  {
    if (fieldTgt->name_ == "U")
    {
      Foam::volVectorField const & fieldSrc =
          mesh_->lookupObject<Foam::volVectorField>(fieldTgt->name_);
      std::vector<double> data(fieldSrc.size() * 3u);
      forAll(mesh_->cells(), cellId)
      {
        data[3 * cellId + 0] = fieldSrc[cellId][0];
        data[3 * cellId + 1] = fieldSrc[cellId][1];
        data[3 * cellId + 2] = fieldSrc[cellId][2];
      }
      fieldTgt->setValues(data, 3u);
    }
    else
    {
      Foam::volScalarField const & fieldSrc =
          mesh_->lookupObject<Foam::volScalarField>(fieldTgt->name_);
      std::vector<double> data(fieldSrc.size());
      forAll(mesh_->cells(), cellId) { data[cellId] = fieldSrc[cellId]; }
      fieldTgt->setValues(data, 1u);
    }
  }
  else if (fieldTgt->mesh_->scope_ == COUPLING_SCOPE::BOUNDARY)
  {
    if (fieldTgt->name_ == "U")
    {
      auto const & fieldSrc =
          mesh_->lookupObject<Foam::volVectorField>(fieldTgt->name_);
      auto const & fieldSrcBd = fieldSrc.boundaryField()[fieldTgt->mesh_->marker_];
      std::vector<double> data(fieldSrcBd.size() * 3u);
      auto count = 0u;
      for (auto const & value: fieldSrcBd)
      {
        data[3 * count + 0] = value[0u];
        data[3 * count + 1] = value[1u];
        data[3 * count + 2] = value[2u];
        count++;
      }
      fieldTgt->setValues(data, 3u);
    }
    else
    {
      auto const & fieldSrc =
          mesh_->lookupObject<Foam::volScalarField>(fieldTgt->name_);
      auto const & fieldSrcBd = fieldSrc.boundaryField()[fieldTgt->mesh_->marker_];
      std::vector<double> data(fieldSrcBd.size());
      auto count = 0u;
      for (auto const & value: fieldSrcBd)
      {
        data[count] = value;
        count++;
      }
      fieldTgt->setValues(data, 1u);
    }
  }
  else
  {
    fmt::println(stderr, "field scope not set");
    std::abort();
  }
}

void ProblemOForg::getFieldData(FieldCoupling const & fieldSrc)
{
  if (fieldSrc.mesh_->scope_ == COUPLING_SCOPE::VOLUME)
  {
    if (fieldSrc.name_ == "U")
    {
      auto & fieldTgt = mesh_->lookupObjectRef<Foam::volVectorField>(fieldSrc.name_);
      for (int i = 0; i < fieldTgt.size(); i++)
      {
        fieldTgt.primitiveFieldRef()[i][0] = fieldSrc[3 * i + 0];
        fieldTgt.primitiveFieldRef()[i][1] = fieldSrc[3 * i + 1];
        fieldTgt.primitiveFieldRef()[i][2] = fieldSrc[3 * i + 2];
      }
    }
    else
    {
      auto & fieldTgt = mesh_->lookupObjectRef<Foam::volScalarField>(fieldSrc.name_);
      for (int i = 0; i < fieldTgt.size(); i++)
        fieldTgt.primitiveFieldRef()[i] = fieldSrc[i];
    }
  }
  else if (fieldSrc.mesh_->scope_ == COUPLING_SCOPE::BOUNDARY)
  {
    if (fieldSrc.name_ == "U")
    {
      auto & fieldTgt = mesh_->lookupObjectRef<Foam::volVectorField>(fieldSrc.name_);
      auto & fieldTgtBd = fieldTgt.boundaryFieldRef()[fieldSrc.mesh_->marker_];
      for (int i = 0; i < fieldTgtBd.size(); i++)
      {
        fieldTgtBd[i][0] = fieldSrc[3 * i + 0];
        fieldTgtBd[i][1] = fieldSrc[3 * i + 1];
        fieldTgtBd[i][2] = fieldSrc[3 * i + 2];
      }
    }
    else
    {
      Foam::volScalarField & fieldTgt =
          mesh_->lookupObjectRef<Foam::volScalarField>(fieldSrc.name_);
      auto & fieldTgtBd = fieldTgt.boundaryFieldRef()[fieldSrc.mesh_->marker_];
      for (int i = 0; i < fieldTgtBd.size(); i++)
        fieldTgtBd[i] = fieldSrc[i];
    }
  }
  else
  {
    fmt::println(stderr, "field scope not set");
    std::abort();
  }
}

bool ProblemOForg::run() const { return pimple_->run(*runTime_); }

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
    // version >12
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

  // TODO: extract the number of iterations
  return 0u;
}

void ProblemOForg::print()
{
  runTime_->write();

  Foam::Info << "ExecutionTime = " << runTime_->elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime_->elapsedClockTime() << " s" << Foam::nl
             << Foam::endl;

  Foam::word const writeStyle = runTime_->controlDict().lookup("writeControl");
  if (writeStyle == "timeStep")
  {
    Foam::label writeInterval =
        runTime_->controlDict().lookupOrDefault("writeInterval", 1);

    if (it % writeInterval == 0)
    {
      printVTK();
    }
  }
  else if (writeStyle == "runTime")
  {
    double const dtWrite =
        runTime_->controlDict().lookupOrDefault("writeInterval", 1.0);
    if (lastPrint_ + dtWrite > time - 1.e-12)
    {
      lastPrint_ = time;
      printVTK();
    }
  }
  else
  {
    fmt::println(stderr, "writeControl {} not supported", writeStyle);
  }
}

void ProblemOForg::printVTK()
{
  for (auto const & [name, field]: fieldsExport_)
  {
    setFieldData(field.get());
    field->printVTK(time, it);
  }
}

Marker ProblemOForg::findRegion(std::string_view name)
{
  auto const patchId = mesh_->boundaryMesh().findPatchID(std::string{name});
  if (patchId < 0)
  {
    fmt::println("region {} not found", name);
    std::abort();
  }

  return patchId;

  // const Foam::fvPatchList & patches = mesh_->boundary();
  // for (auto const & patch: patches)
  // {
  //   if (patch.name() == name)
  //   {
  //     return patch.index();
  //   }
  // }

  // fmt::println("region {} not recognized", name);
  // std::abort();
  // return markerNotSet;
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
