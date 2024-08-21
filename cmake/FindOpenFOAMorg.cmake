# define components
add_library(OpenFOAMorg::OpenFOAM INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::OpenFOAM
  INTERFACE
    $ENV{FOAM_LIBBIN}/libOpenFOAM.so
)
target_include_directories(OpenFOAMorg::OpenFOAM
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude
)

add_library(OpenFOAMorg::MomentumTransportModels_compressible INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::MomentumTransportModels_compressible
  INTERFACE
    $ENV{FOAM_LIBBIN}/libcompressibleMomentumTransportModels.so
)
target_include_directories(OpenFOAMorg::MomentumTransportModels_compressible
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
)

add_library(OpenFOAMorg::MomentumTransportModels_incompressible INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::MomentumTransportModels_incompressible
  INTERFACE
    $ENV{FOAM_LIBBIN}/libincompressibleMomentumTransportModels.so
)
target_include_directories(OpenFOAMorg::MomentumTransportModels_incompressible
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/incompressible/lnInclude
)

add_library(OpenFOAMorg::MomentumTransportModels_momentumTransportModels INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::MomentumTransportModels_momentumTransportModels
  INTERFACE
    $ENV{FOAM_LIBBIN}/libmomentumTransportModels.so
)
target_include_directories(OpenFOAMorg::MomentumTransportModels_momentumTransportModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
)

add_library(OpenFOAMorg::MomentumTransportModels_phaseCompressible INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::MomentumTransportModels_phaseCompressible
  INTERFACE
    $ENV{FOAM_LIBBIN}/libphaseCompressibleMomentumTransportModels.so
)
target_include_directories(OpenFOAMorg::MomentumTransportModels_phaseCompressible
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/phaseCompressible/lnInclude
)

add_library(OpenFOAMorg::MomentumTransportModels_phaseIncompressible INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::MomentumTransportModels_phaseIncompressible
  INTERFACE
    $ENV{FOAM_LIBBIN}/libphaseIncompressibleMomentumTransportModels.so
)
target_include_directories(OpenFOAMorg::MomentumTransportModels_phaseIncompressible
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/phaseIncompressible/lnInclude
)

add_library(OpenFOAMorg::ODE INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::ODE
  INTERFACE
    $ENV{FOAM_LIBBIN}/libODE.so
)
target_include_directories(OpenFOAMorg::ODE
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ODE/lnInclude
)

add_library(OpenFOAMorg::ThermophysicalTransportModels_coupledThermophysicalTransportModels INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_coupledThermophysicalTransportModels
  INTERFACE
    $ENV{FOAM_LIBBIN}/libcoupledThermophysicalTransportModels.so
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_coupledThermophysicalTransportModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/coupledThermophysicalTransportModels/lnInclude
)

add_library(OpenFOAMorg::ThermophysicalTransportModels_fluid INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_fluid
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfluidThermophysicalTransportModel.so
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_fluid
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/fluid/lnInclude
)

add_library(OpenFOAMorg::ThermophysicalTransportModels_fluidMulticomponentThermo INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_fluidMulticomponentThermo
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfluidMulticomponentThermophysicalTransportModels.so
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_fluidMulticomponentThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/fluidMulticomponentThermo/lnInclude
)

add_library(OpenFOAMorg::ThermophysicalTransportModels_fluidThermo INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfluidThermoThermophysicalTransportModels.so
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/fluidThermo/lnInclude
)

add_library(OpenFOAMorg::ThermophysicalTransportModels_phaseFluidMulticomponentThermo INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_phaseFluidMulticomponentThermo
  INTERFACE
    $ENV{FOAM_LIBBIN}/libphaseFluidMulticomponentThermophysicalTransportModels.so
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_phaseFluidMulticomponentThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/phaseFluidMulticomponentThermo/lnInclude
)

add_library(OpenFOAMorg::ThermophysicalTransportModels_phaseFluidThermo INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_phaseFluidThermo
  INTERFACE
    $ENV{FOAM_LIBBIN}/libphaseFluidThermophysicalTransportModels.so
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_phaseFluidThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/phaseFluidThermo/lnInclude
)

add_library(OpenFOAMorg::ThermophysicalTransportModels_solid INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_solid
  INTERFACE
    $ENV{FOAM_LIBBIN}/libsolidThermophysicalTransportModels.so
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_solid
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/solid/lnInclude
)

add_library(OpenFOAMorg::ThermophysicalTransportModels_thermophysicalTransportModel INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_thermophysicalTransportModel
  INTERFACE
    $ENV{FOAM_LIBBIN}/libthermophysicalTransportModel.so
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_thermophysicalTransportModel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
)

add_library(OpenFOAMorg::atmosphericModels INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::atmosphericModels
  INTERFACE
    $ENV{FOAM_LIBBIN}/libatmosphericModels.so
)
target_include_directories(OpenFOAMorg::atmosphericModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/atmosphericModels/lnInclude
)

add_library(OpenFOAMorg::combustionModels INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::combustionModels
  INTERFACE
    $ENV{FOAM_LIBBIN}/libcombustionModels.so
)
target_include_directories(OpenFOAMorg::combustionModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/combustionModels/lnInclude
)

add_library(OpenFOAMorg::conversion INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::conversion
  INTERFACE
    $ENV{FOAM_LIBBIN}/libconversion.so
)
target_include_directories(OpenFOAMorg::conversion
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/conversion/lnInclude
)

add_library(OpenFOAMorg::dynamicMesh INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::dynamicMesh
  INTERFACE
    $ENV{FOAM_LIBBIN}/libdynamicMesh.so
)
target_include_directories(OpenFOAMorg::dynamicMesh
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
)

add_library(OpenFOAMorg::fileFormats INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::fileFormats
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfileFormats.so
)
target_include_directories(OpenFOAMorg::fileFormats
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
)

add_library(OpenFOAMorg::finiteVolume INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::finiteVolume
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfiniteVolume.so
)
target_include_directories(OpenFOAMorg::finiteVolume
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

add_library(OpenFOAMorg::functionObjects_field INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::functionObjects_field
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfieldFunctionObjects.so
)
target_include_directories(OpenFOAMorg::functionObjects_field
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/functionObjects/field/lnInclude
)

add_library(OpenFOAMorg::functionObjects_forces INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::functionObjects_forces
  INTERFACE
    $ENV{FOAM_LIBBIN}/libforces.so
)
target_include_directories(OpenFOAMorg::functionObjects_forces
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/functionObjects/forces/lnInclude
)

add_library(OpenFOAMorg::functionObjects_solvers INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::functionObjects_solvers
  INTERFACE
    $ENV{FOAM_LIBBIN}/libsolverFunctionObjects.so
)
target_include_directories(OpenFOAMorg::functionObjects_solvers
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/functionObjects/solvers/lnInclude
)

add_library(OpenFOAMorg::functionObjects_utilities INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::functionObjects_utilities
  INTERFACE
    $ENV{FOAM_LIBBIN}/libutilityFunctionObjects.so
)
target_include_directories(OpenFOAMorg::functionObjects_utilities
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/functionObjects/utilities/lnInclude
)

add_library(OpenFOAMorg::fvAgglomerationMethods_pairPatchAgglomeration INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::fvAgglomerationMethods_pairPatchAgglomeration
  INTERFACE
    $ENV{FOAM_LIBBIN}/libpairPatchAgglomeration.so
)
target_include_directories(OpenFOAMorg::fvAgglomerationMethods_pairPatchAgglomeration
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fvAgglomerationMethods/pairPatchAgglomeration/lnInclude
)

add_library(OpenFOAMorg::fvConstraints INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::fvConstraints
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfvConstraints.so
)
target_include_directories(OpenFOAMorg::fvConstraints
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fvConstraints/lnInclude
)

add_library(OpenFOAMorg::fvMeshMovers INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::fvMeshMovers
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfvMeshMovers.so
)
target_include_directories(OpenFOAMorg::fvMeshMovers
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fvMeshMovers/lnInclude
)

add_library(OpenFOAMorg::fvMeshStitchers INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::fvMeshStitchers
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfvMeshStitchers.so
)
target_include_directories(OpenFOAMorg::fvMeshStitchers
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fvMeshStitchers/lnInclude
)

add_library(OpenFOAMorg::fvMeshTopoChangers INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::fvMeshTopoChangers
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfvMeshTopoChangers.so
)
target_include_directories(OpenFOAMorg::fvMeshTopoChangers
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fvMeshTopoChangers/lnInclude
)

add_library(OpenFOAMorg::fvMeshTopoChangers_meshToMesh INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::fvMeshTopoChangers_meshToMesh
  INTERFACE
    $ENV{FOAM_LIBBIN}/libmeshToMeshTopoChanger.so
)
target_include_directories(OpenFOAMorg::fvMeshTopoChangers_meshToMesh
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fvMeshTopoChangers/meshToMesh/lnInclude
)

add_library(OpenFOAMorg::fvModels INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::fvModels
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfvModels.so
)
target_include_directories(OpenFOAMorg::fvModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fvModels/lnInclude
)

add_library(OpenFOAMorg::fvMotionSolver INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::fvMotionSolver
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfvMotionSolvers.so
)
target_include_directories(OpenFOAMorg::fvMotionSolver
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fvMotionSolver/lnInclude
)

add_library(OpenFOAMorg::genericPatchFields INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::genericPatchFields
  INTERFACE
    $ENV{FOAM_LIBBIN}/libgenericPatchFields.so
)
target_include_directories(OpenFOAMorg::genericPatchFields
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/genericPatchFields/lnInclude
)

add_library(OpenFOAMorg::genericPatches INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::genericPatches
  INTERFACE
    $ENV{FOAM_LIBBIN}/libgenericPatches.so
)
target_include_directories(OpenFOAMorg::genericPatches
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/genericPatches/lnInclude
)

add_library(OpenFOAMorg::lagrangian_DSMC INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::lagrangian_DSMC
  INTERFACE
    $ENV{FOAM_LIBBIN}/libDSMC.so
)
target_include_directories(OpenFOAMorg::lagrangian_DSMC
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/lagrangian/DSMC/lnInclude
)

add_library(OpenFOAMorg::lagrangian_basic INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::lagrangian_basic
  INTERFACE
    $ENV{FOAM_LIBBIN}/liblagrangian.so
)
target_include_directories(OpenFOAMorg::lagrangian_basic
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/lagrangian/basic/lnInclude
)

add_library(OpenFOAMorg::lagrangian_functionObjects INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::lagrangian_functionObjects
  INTERFACE
    $ENV{FOAM_LIBBIN}/liblagrangianFunctionObjects.so
)
target_include_directories(OpenFOAMorg::lagrangian_functionObjects
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/lagrangian/functionObjects/lnInclude
)

add_library(OpenFOAMorg::lagrangian_molecularDynamics_molecularMeasurements INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::lagrangian_molecularDynamics_molecularMeasurements
  INTERFACE
    $ENV{FOAM_LIBBIN}/libmolecularMeasurements.so
)
target_include_directories(OpenFOAMorg::lagrangian_molecularDynamics_molecularMeasurements
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/lagrangian/molecularDynamics/molecularMeasurements/lnInclude
)

add_library(OpenFOAMorg::lagrangian_molecularDynamics_molecule INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::lagrangian_molecularDynamics_molecule
  INTERFACE
    $ENV{FOAM_LIBBIN}/libmolecule.so
)
target_include_directories(OpenFOAMorg::lagrangian_molecularDynamics_molecule
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/lagrangian/molecularDynamics/molecule/lnInclude
)

add_library(OpenFOAMorg::lagrangian_molecularDynamics_potential INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::lagrangian_molecularDynamics_potential
  INTERFACE
    $ENV{FOAM_LIBBIN}/libpotential.so
)
target_include_directories(OpenFOAMorg::lagrangian_molecularDynamics_potential
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/lagrangian/molecularDynamics/potential/lnInclude
)

add_library(OpenFOAMorg::lagrangian_parcel INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::lagrangian_parcel
  INTERFACE
    $ENV{FOAM_LIBBIN}/liblagrangianParcel.so
)
target_include_directories(OpenFOAMorg::lagrangian_parcel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/lagrangian/parcel/lnInclude
)

add_library(OpenFOAMorg::lagrangian_solidParticle INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::lagrangian_solidParticle
  INTERFACE
    $ENV{FOAM_LIBBIN}/libsolidParticle.so
)
target_include_directories(OpenFOAMorg::lagrangian_solidParticle
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/lagrangian/solidParticle/lnInclude
)

add_library(OpenFOAMorg::mesh_blockMesh INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::mesh_blockMesh
  INTERFACE
    $ENV{FOAM_LIBBIN}/libblockMesh.so
)
target_include_directories(OpenFOAMorg::mesh_blockMesh
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/mesh/blockMesh/lnInclude
)

add_library(OpenFOAMorg::mesh_extrudeModel INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::mesh_extrudeModel
  INTERFACE
    $ENV{FOAM_LIBBIN}/libextrudeModel.so
)
target_include_directories(OpenFOAMorg::mesh_extrudeModel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/mesh/extrudeModel/lnInclude
)

add_library(OpenFOAMorg::mesh_snappyHexMesh INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::mesh_snappyHexMesh
  INTERFACE
    $ENV{FOAM_LIBBIN}/libsnappyHexMesh.so
)
target_include_directories(OpenFOAMorg::mesh_snappyHexMesh
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/mesh/snappyHexMesh/lnInclude
)

add_library(OpenFOAMorg::meshTools INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::meshTools
  INTERFACE
    $ENV{FOAM_LIBBIN}/libmeshTools.so
)
target_include_directories(OpenFOAMorg::meshTools
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

add_library(OpenFOAMorg::multiphaseModels_multiphaseProperties INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::multiphaseModels_multiphaseProperties
  INTERFACE
    $ENV{FOAM_LIBBIN}/libmultiphaseProperties.so
)
target_include_directories(OpenFOAMorg::multiphaseModels_multiphaseProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/multiphaseModels/multiphaseProperties/lnInclude
)

add_library(OpenFOAMorg::parallel_decompose_decompositionMethods INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::parallel_decompose_decompositionMethods
  INTERFACE
    $ENV{FOAM_LIBBIN}/libdecompositionMethods.so
)
target_include_directories(OpenFOAMorg::parallel_decompose_decompositionMethods
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/parallel/decompose/decompositionMethods/lnInclude
)

add_library(OpenFOAMorg::parallel_distributed INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::parallel_distributed
  INTERFACE
    $ENV{FOAM_LIBBIN}/libdistributed.so
)
target_include_directories(OpenFOAMorg::parallel_distributed
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/parallel/distributed/lnInclude
)

add_library(OpenFOAMorg::parallel_parallel INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::parallel_parallel
  INTERFACE
    $ENV{FOAM_LIBBIN}/libparallel.so
)
target_include_directories(OpenFOAMorg::parallel_parallel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/parallel/parallel/lnInclude
)

add_library(OpenFOAMorg::physicalProperties INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::physicalProperties
  INTERFACE
    $ENV{FOAM_LIBBIN}/libphysicalProperties.so
)
target_include_directories(OpenFOAMorg::physicalProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
)

add_library(OpenFOAMorg::radiationModels INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::radiationModels
  INTERFACE
    $ENV{FOAM_LIBBIN}/libradiationModels.so
)
target_include_directories(OpenFOAMorg::radiationModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/radiationModels/lnInclude
)

add_library(OpenFOAMorg::randomProcesses INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::randomProcesses
  INTERFACE
    $ENV{FOAM_LIBBIN}/librandomProcesses.so
)
target_include_directories(OpenFOAMorg::randomProcesses
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/randomProcesses/lnInclude
)

add_library(OpenFOAMorg::renumber_renumberMethods INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::renumber_renumberMethods
  INTERFACE
    $ENV{FOAM_LIBBIN}/librenumberMethods.so
)
target_include_directories(OpenFOAMorg::renumber_renumberMethods
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/renumber/renumberMethods/lnInclude
)

add_library(OpenFOAMorg::rigidBodyDynamics INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::rigidBodyDynamics
  INTERFACE
    $ENV{FOAM_LIBBIN}/librigidBodyDynamics.so
)
target_include_directories(OpenFOAMorg::rigidBodyDynamics
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/rigidBodyDynamics/lnInclude
)

add_library(OpenFOAMorg::rigidBodyMeshMotion INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::rigidBodyMeshMotion
  INTERFACE
    $ENV{FOAM_LIBBIN}/librigidBodyMeshMotion.so
)
target_include_directories(OpenFOAMorg::rigidBodyMeshMotion
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/rigidBodyMeshMotion/lnInclude
)

add_library(OpenFOAMorg::rigidBodyState INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::rigidBodyState
  INTERFACE
    $ENV{FOAM_LIBBIN}/librigidBodyState.so
)
target_include_directories(OpenFOAMorg::rigidBodyState
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/rigidBodyState/lnInclude
)

add_library(OpenFOAMorg::sampling INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::sampling
  INTERFACE
    $ENV{FOAM_LIBBIN}/libsampling.so
)
target_include_directories(OpenFOAMorg::sampling
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/sampling/lnInclude
)

add_library(OpenFOAMorg::sixDoFRigidBodyMotion INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::sixDoFRigidBodyMotion
  INTERFACE
    $ENV{FOAM_LIBBIN}/libsixDoFRigidBodyMotion.so
)
target_include_directories(OpenFOAMorg::sixDoFRigidBodyMotion
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/sixDoFRigidBodyMotion/lnInclude
)

add_library(OpenFOAMorg::sixDoFRigidBodyState INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::sixDoFRigidBodyState
  INTERFACE
    $ENV{FOAM_LIBBIN}/libsixDoFRigidBodyState.so
)
target_include_directories(OpenFOAMorg::sixDoFRigidBodyState
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/sixDoFRigidBodyState/lnInclude
)

add_library(OpenFOAMorg::specieTransfer INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::specieTransfer
  INTERFACE
    $ENV{FOAM_LIBBIN}/libspecieTransfer.so
)
target_include_directories(OpenFOAMorg::specieTransfer
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/specieTransfer/lnInclude
)

add_library(OpenFOAMorg::surfMesh INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::surfMesh
  INTERFACE
    $ENV{FOAM_LIBBIN}/libsurfMesh.so
)
target_include_directories(OpenFOAMorg::surfMesh
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/surfMesh/lnInclude
)

add_library(OpenFOAMorg::thermophysicalModels_barotropicCompressibilityModel INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::thermophysicalModels_barotropicCompressibilityModel
  INTERFACE
    $ENV{FOAM_LIBBIN}/libbarotropicCompressibilityModel.so
)
target_include_directories(OpenFOAMorg::thermophysicalModels_barotropicCompressibilityModel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/barotropicCompressibilityModel/lnInclude
)

add_library(OpenFOAMorg::thermophysicalModels_basic INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::thermophysicalModels_basic
  INTERFACE
    $ENV{FOAM_LIBBIN}/libfluidThermophysicalModels.so
)
target_include_directories(OpenFOAMorg::thermophysicalModels_basic
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
)

add_library(OpenFOAMorg::thermophysicalModels_chemistryModel INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::thermophysicalModels_chemistryModel
  INTERFACE
    $ENV{FOAM_LIBBIN}/libchemistryModel.so
)
target_include_directories(OpenFOAMorg::thermophysicalModels_chemistryModel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/chemistryModel/lnInclude
)

add_library(OpenFOAMorg::thermophysicalModels_ignition INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::thermophysicalModels_ignition
  INTERFACE
    $ENV{FOAM_LIBBIN}/libXiIgnition.so
)
target_include_directories(OpenFOAMorg::thermophysicalModels_ignition
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/ignition/lnInclude
)

add_library(OpenFOAMorg::thermophysicalModels_laminarFlameSpeed INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::thermophysicalModels_laminarFlameSpeed
  INTERFACE
    $ENV{FOAM_LIBBIN}/liblaminarFlameSpeedModels.so
)
target_include_directories(OpenFOAMorg::thermophysicalModels_laminarFlameSpeed
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/laminarFlameSpeed/lnInclude
)

add_library(OpenFOAMorg::thermophysicalModels_multicomponentThermo INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::thermophysicalModels_multicomponentThermo
  INTERFACE
    $ENV{FOAM_LIBBIN}/libmulticomponentThermophysicalModels.so
)
target_include_directories(OpenFOAMorg::thermophysicalModels_multicomponentThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/multicomponentThermo/lnInclude
)

add_library(OpenFOAMorg::thermophysicalModels_saturationModels INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::thermophysicalModels_saturationModels
  INTERFACE
    $ENV{FOAM_LIBBIN}/libsaturationModels.so
)
target_include_directories(OpenFOAMorg::thermophysicalModels_saturationModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/saturationModels/lnInclude
)

add_library(OpenFOAMorg::thermophysicalModels_solidThermo INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::thermophysicalModels_solidThermo
  INTERFACE
    $ENV{FOAM_LIBBIN}/libsolidThermo.so
)
target_include_directories(OpenFOAMorg::thermophysicalModels_solidThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/solidThermo/lnInclude
)

add_library(OpenFOAMorg::thermophysicalModels_specie INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::thermophysicalModels_specie
  INTERFACE
    $ENV{FOAM_LIBBIN}/libspecie.so
)
target_include_directories(OpenFOAMorg::thermophysicalModels_specie
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
)

add_library(OpenFOAMorg::thermophysicalModels_thermophysicalProperties INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::thermophysicalModels_thermophysicalProperties
  INTERFACE
    $ENV{FOAM_LIBBIN}/libthermophysicalProperties.so
)
target_include_directories(OpenFOAMorg::thermophysicalModels_thermophysicalProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/thermophysicalProperties/lnInclude
)

add_library(OpenFOAMorg::triSurface INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::triSurface
  INTERFACE
    $ENV{FOAM_LIBBIN}/libtriSurface.so
)
target_include_directories(OpenFOAMorg::triSurface
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/triSurface/lnInclude
)

add_library(OpenFOAMorg::twoPhaseModels_compressibleCavitation INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::twoPhaseModels_compressibleCavitation
  INTERFACE
    $ENV{FOAM_LIBBIN}/libcompressibleCavitationModels.so
)
target_include_directories(OpenFOAMorg::twoPhaseModels_compressibleCavitation
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/compressibleCavitation/lnInclude
)

add_library(OpenFOAMorg::twoPhaseModels_compressibleInterfaceProperties INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::twoPhaseModels_compressibleInterfaceProperties
  INTERFACE
    $ENV{FOAM_LIBBIN}/libcompressibleInterfaceProperties.so
)
target_include_directories(OpenFOAMorg::twoPhaseModels_compressibleInterfaceProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/compressibleInterfaceProperties/lnInclude
)

add_library(OpenFOAMorg::twoPhaseModels_compressibleTwoPhases INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::twoPhaseModels_compressibleTwoPhases
  INTERFACE
    $ENV{FOAM_LIBBIN}/libcompressibleTwoPhases.so
)
target_include_directories(OpenFOAMorg::twoPhaseModels_compressibleTwoPhases
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/compressibleTwoPhases/lnInclude
)

add_library(OpenFOAMorg::twoPhaseModels_incompressibleCavitation INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::twoPhaseModels_incompressibleCavitation
  INTERFACE
    $ENV{FOAM_LIBBIN}/libincompressibleCavitationModels.so
)
target_include_directories(OpenFOAMorg::twoPhaseModels_incompressibleCavitation
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/incompressibleCavitation/lnInclude
)

add_library(OpenFOAMorg::twoPhaseModels_incompressibleTwoPhases INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::twoPhaseModels_incompressibleTwoPhases
  INTERFACE
    $ENV{FOAM_LIBBIN}/libincompressibleTwoPhases.so
)
target_include_directories(OpenFOAMorg::twoPhaseModels_incompressibleTwoPhases
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/incompressibleTwoPhases/lnInclude
)

add_library(OpenFOAMorg::twoPhaseModels_interfaceCompression INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::twoPhaseModels_interfaceCompression
  INTERFACE
    $ENV{FOAM_LIBBIN}/libinterfaceCompression.so
)
target_include_directories(OpenFOAMorg::twoPhaseModels_interfaceCompression
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/interfaceCompression/lnInclude
)

add_library(OpenFOAMorg::twoPhaseModels_interfaceProperties INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::twoPhaseModels_interfaceProperties
  INTERFACE
    $ENV{FOAM_LIBBIN}/libinterfaceProperties.so
)
target_include_directories(OpenFOAMorg::twoPhaseModels_interfaceProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/interfaceProperties/lnInclude
)

add_library(OpenFOAMorg::twoPhaseModels_twoPhaseMixture INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::twoPhaseModels_twoPhaseMixture
  INTERFACE
    $ENV{FOAM_LIBBIN}/libtwoPhaseMixture.so
)
target_include_directories(OpenFOAMorg::twoPhaseModels_twoPhaseMixture
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/twoPhaseMixture/lnInclude
)

add_library(OpenFOAMorg::twoPhaseModels_twoPhaseProperties INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::twoPhaseModels_twoPhaseProperties
  INTERFACE
    $ENV{FOAM_LIBBIN}/libtwoPhaseProperties.so
)
target_include_directories(OpenFOAMorg::twoPhaseModels_twoPhaseProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/twoPhaseProperties/lnInclude
)

add_library(OpenFOAMorg::waves INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::waves
  INTERFACE
    $ENV{FOAM_LIBBIN}/libwaves.so
)
target_include_directories(OpenFOAMorg::waves
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/waves/lnInclude
)

# set up dependencies
target_link_libraries(OpenFOAMorg::OpenFOAM
  INTERFACE
  $ENV{WM_PROJECT_DIR}/platforms/linux64GccDPInt32Opt/lib/openmpi-system/libPstream.so
)
target_include_directories(OpenFOAMorg::OpenFOAM
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/OSspecific/POSIX/lnInclude
)

target_link_libraries(OpenFOAMorg::MomentumTransportModels_compressible
  INTERFACE
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::MomentumTransportModels_compressible
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::MomentumTransportModels_incompressible
  INTERFACE
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::MomentumTransportModels_incompressible
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::MomentumTransportModels_momentumTransportModels
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::MomentumTransportModels_momentumTransportModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::MomentumTransportModels_phaseCompressible
  INTERFACE
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::MomentumTransportModels_phaseCompressible
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::MomentumTransportModels_phaseIncompressible
  INTERFACE
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::MomentumTransportModels_phaseIncompressible
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/incompressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

# component ODE has no dependencies

target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_coupledThermophysicalTransportModels
  INTERFACE
    OpenFOAMorg::ThermophysicalTransportModels_thermophysicalTransportModel
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_coupledThermophysicalTransportModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_fluid
  INTERFACE
    OpenFOAMorg::ThermophysicalTransportModels_thermophysicalTransportModel
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_solidThermo
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_fluid
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/solidThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_fluidMulticomponentThermo
  INTERFACE
    OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_multicomponentThermo
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_fluidMulticomponentThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/fluid/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/multicomponentThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
  INTERFACE
    OpenFOAMorg::ThermophysicalTransportModels_fluid
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_solidThermo
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/fluid/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/solidThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_phaseFluidMulticomponentThermo
  INTERFACE
    OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_multicomponentThermo
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_phaseFluidMulticomponentThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/fluid/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/phaseCompressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/multicomponentThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_phaseFluidThermo
  INTERFACE
    OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_phaseFluidThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/fluid/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/phaseCompressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_solid
  INTERFACE
    OpenFOAMorg::ThermophysicalTransportModels_thermophysicalTransportModel
    OpenFOAMorg::thermophysicalModels_solidThermo
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_solid
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/solidThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::ThermophysicalTransportModels_thermophysicalTransportModel
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::ThermophysicalTransportModels_thermophysicalTransportModel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::atmosphericModels
  INTERFACE
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::MomentumTransportModels_incompressible
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
    OpenFOAMorg::sampling
    OpenFOAMorg::triSurface
    OpenFOAMorg::fvModels
    OpenFOAMorg::fvConstraints
)
target_include_directories(OpenFOAMorg::atmosphericModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/incompressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/sampling/lnInclude
    $ENV{WM_PROJECT_DIR}/src/triSurface/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fvModels/lnInclude
)

target_link_libraries(OpenFOAMorg::combustionModels
  INTERFACE
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::MomentumTransportModels_compressible
    OpenFOAMorg::thermophysicalModels_chemistryModel
    OpenFOAMorg::radiationModels
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::combustionModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/multicomponentThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/chemistryModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/radiationModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::conversion
  INTERFACE
    OpenFOAMorg::fileFormats
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::conversion
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::dynamicMesh
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::triSurface
    OpenFOAMorg::mesh_extrudeModel
)
target_include_directories(OpenFOAMorg::dynamicMesh
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/triSurface/lnInclude
    $ENV{WM_PROJECT_DIR}/src/mesh/extrudeModel/lnInclude
)

# component fileFormats has no dependencies

target_link_libraries(OpenFOAMorg::finiteVolume
  INTERFACE
    OpenFOAMorg::OpenFOAM
    OpenFOAMorg::triSurface
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::finiteVolume
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/triSurface/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::functionObjects_field
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_solidThermo
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::MomentumTransportModels_incompressible
    OpenFOAMorg::MomentumTransportModels_compressible
    OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
    OpenFOAMorg::ThermophysicalTransportModels_solid
    OpenFOAMorg::meshTools
    OpenFOAMorg::surfMesh
    OpenFOAMorg::lagrangian_basic
    OpenFOAMorg::fileFormats
    OpenFOAMorg::sampling
    OpenFOAMorg::surfMesh
)
target_include_directories(OpenFOAMorg::functionObjects_field
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/lagrangian/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
    $ENV{WM_PROJECT_DIR}/src/sampling/lnInclude
    $ENV{WM_PROJECT_DIR}/src/surfMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/solidThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/incompressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/fluid/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/solid/lnInclude
)

target_link_libraries(OpenFOAMorg::functionObjects_forces
  INTERFACE
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::MomentumTransportModels_incompressible
    OpenFOAMorg::MomentumTransportModels_compressible
    OpenFOAMorg::MomentumTransportModels_phaseIncompressible
    OpenFOAMorg::MomentumTransportModels_phaseCompressible
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::fileFormats
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::functionObjects_forces
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/incompressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/phaseIncompressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/phaseCompressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::functionObjects_solvers
  INTERFACE
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::twoPhaseModels_interfaceCompression
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::MomentumTransportModels_incompressible
    OpenFOAMorg::MomentumTransportModels_compressible
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
    OpenFOAMorg::fvModels
    OpenFOAMorg::fvConstraints
)
target_include_directories(OpenFOAMorg::functionObjects_solvers
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/interfaceCompression/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/incompressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::functionObjects_utilities
  INTERFACE
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::functionObjects_utilities
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::fvAgglomerationMethods_pairPatchAgglomeration
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::OpenFOAM
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::fvAgglomerationMethods_pairPatchAgglomeration
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude
)

target_link_libraries(OpenFOAMorg::fvConstraints
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::sampling
    OpenFOAMorg::meshTools
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::MomentumTransportModels_compressible
    OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
)
target_include_directories(OpenFOAMorg::fvConstraints
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/sampling/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/solidThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/fluid/lnInclude
)

target_link_libraries(OpenFOAMorg::fvMeshMovers
  INTERFACE
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
    OpenFOAMorg::fvMeshStitchers
)
target_include_directories(OpenFOAMorg::fvMeshMovers
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::fvMeshStitchers
  INTERFACE
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::fvMeshStitchers
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::fvMeshTopoChangers
  INTERFACE
    OpenFOAMorg::triSurface
    OpenFOAMorg::meshTools
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::fvMeshStitchers
)
target_include_directories(OpenFOAMorg::fvMeshTopoChangers
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/triSurface/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::fvMeshTopoChangers_meshToMesh
  INTERFACE
    OpenFOAMorg::triSurface
    OpenFOAMorg::meshTools
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::fvMeshTopoChangers_meshToMesh
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/triSurface/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::fvModels
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::sampling
    OpenFOAMorg::meshTools
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::MomentumTransportModels_compressible
    OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
)
target_include_directories(OpenFOAMorg::fvModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/sampling/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/solidThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/fluid/lnInclude
)

target_link_libraries(OpenFOAMorg::fvMotionSolver
  INTERFACE
    OpenFOAMorg::triSurface
    OpenFOAMorg::meshTools
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::fvMeshMovers
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::fileFormats
)
target_include_directories(OpenFOAMorg::fvMotionSolver
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/triSurface/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fvMeshMovers/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
    $ENV{WM_PROJECT_DIR}/src/functionObjects/forces/lnInclude
)

target_link_libraries(OpenFOAMorg::genericPatchFields
  INTERFACE
    OpenFOAMorg::genericPatches
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::genericPatchFields
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/genericPatches/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

# component genericPatches has no dependencies

target_link_libraries(OpenFOAMorg::lagrangian_DSMC
  INTERFACE
    OpenFOAMorg::lagrangian_basic
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::lagrangian_DSMC
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/lagrangian/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::lagrangian_basic
  INTERFACE
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::lagrangian_basic
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::lagrangian_functionObjects
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::meshTools
    OpenFOAMorg::lagrangian_basic
    OpenFOAMorg::lagrangian_parcel
    OpenFOAMorg::functionObjects_utilities
)
target_include_directories(OpenFOAMorg::lagrangian_functionObjects
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/lagrangian/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/lagrangian/parcel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/lagrangian/DSMC/lnInclude
    $ENV{WM_PROJECT_DIR}/src/functionObjects/utilities/lnInclude
)

# component lagrangian_molecularDynamics_molecularMeasurements has no dependencies

target_link_libraries(OpenFOAMorg::lagrangian_molecularDynamics_molecule
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
    OpenFOAMorg::lagrangian_basic
    OpenFOAMorg::lagrangian_molecularDynamics_potential
    OpenFOAMorg::lagrangian_molecularDynamics_molecularMeasurements
)
target_include_directories(OpenFOAMorg::lagrangian_molecularDynamics_molecule
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/lagrangian/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/lagrangian/molecularDynamics/potential/lnInclude
    $ENV{WM_PROJECT_DIR}/src/lagrangian/molecularDynamics/molecularMeasurements/lnInclude
)

target_link_libraries(OpenFOAMorg::lagrangian_molecularDynamics_potential
  INTERFACE
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::lagrangian_molecularDynamics_potential
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/lagrangian
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::lagrangian_parcel
  INTERFACE
    OpenFOAMorg::lagrangian_basic
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_thermophysicalProperties
    OpenFOAMorg::thermophysicalModels_multicomponentThermo
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::MomentumTransportModels_incompressible
    OpenFOAMorg::MomentumTransportModels_compressible
    OpenFOAMorg::radiationModels
    OpenFOAMorg::sampling
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::lagrangian_parcel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/lagrangian/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/thermophysicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/multicomponentThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/incompressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/radiationModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/sampling/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::lagrangian_solidParticle
  INTERFACE
    OpenFOAMorg::lagrangian_basic
    OpenFOAMorg::meshTools
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::lagrangian_solidParticle
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/lagrangian/basic/lnInclude
)

target_link_libraries(OpenFOAMorg::mesh_blockMesh
  INTERFACE
    OpenFOAMorg::meshTools
    OpenFOAMorg::fileFormats
    OpenFOAMorg::surfMesh
)
target_include_directories(OpenFOAMorg::mesh_blockMesh
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
    $ENV{WM_PROJECT_DIR}/src/surfMesh/lnInclude
)

target_link_libraries(OpenFOAMorg::mesh_extrudeModel
  INTERFACE
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::mesh_extrudeModel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
)

target_link_libraries(OpenFOAMorg::mesh_snappyHexMesh
  INTERFACE
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::lagrangian_basic
    OpenFOAMorg::meshTools
    OpenFOAMorg::fileFormats
    OpenFOAMorg::surfMesh
    OpenFOAMorg::triSurface
    OpenFOAMorg::parallel_distributed
)
target_include_directories(OpenFOAMorg::mesh_snappyHexMesh
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/parallel/decompose/decompositionMethods/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/lagrangian/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
    $ENV{WM_PROJECT_DIR}/src/surfMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/triSurface/lnInclude
)

target_link_libraries(OpenFOAMorg::meshTools
  INTERFACE
    OpenFOAMorg::triSurface
    OpenFOAMorg::surfMesh
    OpenFOAMorg::fileFormats
)
target_include_directories(OpenFOAMorg::meshTools
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/triSurface/lnInclude
    $ENV{WM_PROJECT_DIR}/src/surfMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
)

target_link_libraries(OpenFOAMorg::multiphaseModels_multiphaseProperties
  INTERFACE
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::multiphaseModels_multiphaseProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::parallel_decompose_decompositionMethods
  INTERFACE
    OpenFOAMorg::meshTools
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::parallel_decompose_decompositionMethods
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::parallel_distributed
  INTERFACE
    OpenFOAMorg::triSurface
    OpenFOAMorg::parallel_decompose_decompositionMethods
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::parallel_distributed
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/triSurface/lnInclude
    $ENV{WM_PROJECT_DIR}/src/parallel/decompose/decompositionMethods/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::parallel_parallel
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
    OpenFOAMorg::parallel_decompose_decompositionMethods
    OpenFOAMorg::dynamicMesh
)
target_include_directories(OpenFOAMorg::parallel_parallel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/parallel/decompose/decompositionMethods/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
)

target_link_libraries(OpenFOAMorg::physicalProperties
  INTERFACE
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::physicalProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::radiationModels
  INTERFACE
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::thermophysicalModels_solidThermo
    OpenFOAMorg::thermophysicalModels_thermophysicalProperties
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::radiationModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/solidThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/thermophysicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/multicomponentThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::randomProcesses
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::sampling
)
target_include_directories(OpenFOAMorg::randomProcesses
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/sampling/lnInclude
)

target_link_libraries(OpenFOAMorg::renumber_renumberMethods
  INTERFACE
    OpenFOAMorg::parallel_decompose_decompositionMethods
    OpenFOAMorg::dynamicMesh
)
target_include_directories(OpenFOAMorg::renumber_renumberMethods
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/parallel/decompose/decompositionMethods/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

# component rigidBodyDynamics has no dependencies

target_link_libraries(OpenFOAMorg::rigidBodyMeshMotion
  INTERFACE
    OpenFOAMorg::rigidBodyDynamics
    OpenFOAMorg::functionObjects_forces
    OpenFOAMorg::meshTools
    OpenFOAMorg::fileFormats
    OpenFOAMorg::dynamicMesh
)
target_include_directories(OpenFOAMorg::rigidBodyMeshMotion
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/rigidBodyDynamics/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/functionObjects/forces/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
)

target_link_libraries(OpenFOAMorg::rigidBodyState
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::fvMeshMovers
    OpenFOAMorg::rigidBodyDynamics
)
target_include_directories(OpenFOAMorg::rigidBodyState
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fvMeshMovers/lnInclude
    $ENV{WM_PROJECT_DIR}/src/rigidBodyDynamics/lnInclude
)

target_link_libraries(OpenFOAMorg::sampling
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
    OpenFOAMorg::surfMesh
    OpenFOAMorg::fileFormats
    OpenFOAMorg::triSurface
    OpenFOAMorg::lagrangian_basic
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::conversion
)
target_include_directories(OpenFOAMorg::sampling
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/surfMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
    $ENV{WM_PROJECT_DIR}/src/triSurface/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/conversion/lnInclude
    $ENV{WM_PROJECT_DIR}/src/lagrangian/basic/lnInclude
)

target_link_libraries(OpenFOAMorg::sixDoFRigidBodyMotion
  INTERFACE
    OpenFOAMorg::functionObjects_forces
    OpenFOAMorg::meshTools
    OpenFOAMorg::fileFormats
    OpenFOAMorg::dynamicMesh
)
target_include_directories(OpenFOAMorg::sixDoFRigidBodyMotion
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/functionObjects/forces/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
)

target_link_libraries(OpenFOAMorg::sixDoFRigidBodyState
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::fvMeshMovers
    OpenFOAMorg::sixDoFRigidBodyMotion
)
target_include_directories(OpenFOAMorg::sixDoFRigidBodyState
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/fvMeshMovers/lnInclude
    $ENV{WM_PROJECT_DIR}/src/sixDoFRigidBodyMotion/lnInclude
)

target_link_libraries(OpenFOAMorg::specieTransfer
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::thermophysicalModels_multicomponentThermo
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::MomentumTransportModels_compressible
    OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
)
target_include_directories(OpenFOAMorg::specieTransfer
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/multicomponentThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/thermophysicalTransportModel/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/fluid/lnInclude
)

target_link_libraries(OpenFOAMorg::surfMesh
  INTERFACE
    OpenFOAMorg::fileFormats
)
target_include_directories(OpenFOAMorg::surfMesh
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
)

target_link_libraries(OpenFOAMorg::thermophysicalModels_barotropicCompressibilityModel
  INTERFACE
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::thermophysicalModels_barotropicCompressibilityModel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::thermophysicalModels_basic
  INTERFACE
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::thermophysicalModels_thermophysicalProperties
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::thermophysicalModels_basic
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/thermophysicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::thermophysicalModels_chemistryModel
  INTERFACE
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_multicomponentThermo
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::ODE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::thermophysicalModels_chemistryModel
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/multicomponentThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/ODE/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::thermophysicalModels_ignition
  INTERFACE
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::thermophysicalModels_ignition
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::thermophysicalModels_laminarFlameSpeed
  INTERFACE
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::thermophysicalModels_laminarFlameSpeed
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/multicomponentThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::thermophysicalModels_multicomponentThermo
  INTERFACE
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::thermophysicalModels_multicomponentThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::thermophysicalModels_saturationModels
  INTERFACE
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_multicomponentThermo
    OpenFOAMorg::thermophysicalModels_specie
)
target_include_directories(OpenFOAMorg::thermophysicalModels_saturationModels
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/thermophysicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/multicomponentThermo/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/phaseCompressible/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/sampling/lnInclude
)

target_link_libraries(OpenFOAMorg::thermophysicalModels_solidThermo
  INTERFACE
    OpenFOAMorg::meshTools
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_multicomponentThermo
)
target_include_directories(OpenFOAMorg::thermophysicalModels_solidThermo
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/multicomponentThermo/lnInclude
)

# component thermophysicalModels_specie has no dependencies

target_link_libraries(OpenFOAMorg::thermophysicalModels_thermophysicalProperties
  INTERFACE
    OpenFOAMorg::thermophysicalModels_specie
)
target_include_directories(OpenFOAMorg::thermophysicalModels_thermophysicalProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
)

target_link_libraries(OpenFOAMorg::triSurface
  INTERFACE
    OpenFOAMorg::fileFormats
    OpenFOAMorg::surfMesh
)
target_include_directories(OpenFOAMorg::triSurface
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/fileFormats/lnInclude
    $ENV{WM_PROJECT_DIR}/src/surfMesh/lnInclude
)

target_link_libraries(OpenFOAMorg::twoPhaseModels_compressibleCavitation
  INTERFACE
    OpenFOAMorg::twoPhaseModels_compressibleTwoPhases
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_saturationModels
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::twoPhaseModels_compressibleCavitation
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/twoPhaseMixture/lnInclude
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/compressibleTwoPhases/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/saturationModels/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::twoPhaseModels_compressibleInterfaceProperties
  INTERFACE
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::thermophysicalModels_thermophysicalProperties
    OpenFOAMorg::twoPhaseModels_twoPhaseProperties
    OpenFOAMorg::twoPhaseModels_interfaceProperties
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::twoPhaseModels_compressibleInterfaceProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/thermophysicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/twoPhaseProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/interfaceProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::twoPhaseModels_compressibleTwoPhases
  INTERFACE
    OpenFOAMorg::twoPhaseModels_twoPhaseMixture
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::twoPhaseModels_compressibleTwoPhases
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/twoPhaseMixture/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::twoPhaseModels_incompressibleCavitation
  INTERFACE
    OpenFOAMorg::twoPhaseModels_incompressibleTwoPhases
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
)
target_include_directories(OpenFOAMorg::twoPhaseModels_incompressibleCavitation
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/twoPhaseMixture/lnInclude
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/incompressibleTwoPhases/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
)

target_link_libraries(OpenFOAMorg::twoPhaseModels_incompressibleTwoPhases
  INTERFACE
    OpenFOAMorg::twoPhaseModels_twoPhaseMixture
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::twoPhaseModels_incompressibleTwoPhases
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/twoPhaseModels/twoPhaseMixture/lnInclude
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::twoPhaseModels_interfaceCompression
  INTERFACE
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::twoPhaseModels_interfaceCompression
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::twoPhaseModels_interfaceProperties
  INTERFACE
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::twoPhaseModels_interfaceProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::twoPhaseModels_twoPhaseMixture
  INTERFACE
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::twoPhaseModels_twoPhaseMixture
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/physicalProperties/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::twoPhaseModels_twoPhaseProperties
  INTERFACE
    OpenFOAMorg::finiteVolume
)
target_include_directories(OpenFOAMorg::twoPhaseModels_twoPhaseProperties
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
)

target_link_libraries(OpenFOAMorg::waves
  INTERFACE
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::meshTools
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::atmosphericModels
)
target_include_directories(OpenFOAMorg::waves
  INTERFACE
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude
    $ENV{WM_PROJECT_DIR}/src/dynamicMesh/lnInclude
    $ENV{WM_PROJECT_DIR}/src/atmosphericModels/lnInclude
)

# definitions
target_compile_definitions(OpenFOAMorg::OpenFOAM
  INTERFACE
    LIB_NAME=libNULL.so
    linux64
    WM_ARCH_OPTION=64
    WM_DP
    WM_LABEL_SIZE=32
    NoRepository
)

# add ::ALL target
add_library(OpenFOAMorg::ALL INTERFACE IMPORTED)
target_link_libraries(OpenFOAMorg::ALL
  INTERFACE
    OpenFOAMorg::OpenFOAM
    OpenFOAMorg::MomentumTransportModels_compressible
    OpenFOAMorg::MomentumTransportModels_incompressible
    OpenFOAMorg::MomentumTransportModels_momentumTransportModels
    OpenFOAMorg::MomentumTransportModels_phaseCompressible
    OpenFOAMorg::MomentumTransportModels_phaseIncompressible
    OpenFOAMorg::ODE
    OpenFOAMorg::ThermophysicalTransportModels_coupledThermophysicalTransportModels
    OpenFOAMorg::ThermophysicalTransportModels_fluid
    OpenFOAMorg::ThermophysicalTransportModels_fluidMulticomponentThermo
    OpenFOAMorg::ThermophysicalTransportModels_fluidThermo
    OpenFOAMorg::ThermophysicalTransportModels_phaseFluidMulticomponentThermo
    OpenFOAMorg::ThermophysicalTransportModels_phaseFluidThermo
    OpenFOAMorg::ThermophysicalTransportModels_solid
    OpenFOAMorg::ThermophysicalTransportModels_thermophysicalTransportModel
    OpenFOAMorg::atmosphericModels
    OpenFOAMorg::combustionModels
    OpenFOAMorg::conversion
    OpenFOAMorg::dynamicMesh
    OpenFOAMorg::fileFormats
    OpenFOAMorg::finiteVolume
    OpenFOAMorg::functionObjects_field
    OpenFOAMorg::functionObjects_forces
    OpenFOAMorg::functionObjects_solvers
    OpenFOAMorg::functionObjects_utilities
    OpenFOAMorg::fvAgglomerationMethods_pairPatchAgglomeration
    OpenFOAMorg::fvConstraints
    OpenFOAMorg::fvMeshMovers
    OpenFOAMorg::fvMeshStitchers
    OpenFOAMorg::fvMeshTopoChangers
    OpenFOAMorg::fvMeshTopoChangers_meshToMesh
    OpenFOAMorg::fvModels
    OpenFOAMorg::fvMotionSolver
    OpenFOAMorg::genericPatchFields
    OpenFOAMorg::genericPatches
    OpenFOAMorg::lagrangian_DSMC
    OpenFOAMorg::lagrangian_basic
    OpenFOAMorg::lagrangian_functionObjects
    OpenFOAMorg::lagrangian_molecularDynamics_molecularMeasurements
    OpenFOAMorg::lagrangian_molecularDynamics_molecule
    OpenFOAMorg::lagrangian_molecularDynamics_potential
    OpenFOAMorg::lagrangian_parcel
    OpenFOAMorg::lagrangian_solidParticle
    OpenFOAMorg::mesh_blockMesh
    OpenFOAMorg::mesh_extrudeModel
    OpenFOAMorg::mesh_snappyHexMesh
    OpenFOAMorg::meshTools
    OpenFOAMorg::multiphaseModels_multiphaseProperties
    OpenFOAMorg::parallel_decompose_decompositionMethods
    OpenFOAMorg::parallel_distributed
    OpenFOAMorg::parallel_parallel
    OpenFOAMorg::physicalProperties
    OpenFOAMorg::radiationModels
    OpenFOAMorg::randomProcesses
    OpenFOAMorg::renumber_renumberMethods
    OpenFOAMorg::rigidBodyDynamics
    OpenFOAMorg::rigidBodyMeshMotion
    OpenFOAMorg::rigidBodyState
    OpenFOAMorg::sampling
    OpenFOAMorg::sixDoFRigidBodyMotion
    OpenFOAMorg::sixDoFRigidBodyState
    OpenFOAMorg::specieTransfer
    OpenFOAMorg::surfMesh
    OpenFOAMorg::thermophysicalModels_barotropicCompressibilityModel
    OpenFOAMorg::thermophysicalModels_basic
    OpenFOAMorg::thermophysicalModels_chemistryModel
    OpenFOAMorg::thermophysicalModels_ignition
    OpenFOAMorg::thermophysicalModels_laminarFlameSpeed
    OpenFOAMorg::thermophysicalModels_multicomponentThermo
    OpenFOAMorg::thermophysicalModels_saturationModels
    OpenFOAMorg::thermophysicalModels_solidThermo
    OpenFOAMorg::thermophysicalModels_specie
    OpenFOAMorg::thermophysicalModels_thermophysicalProperties
    OpenFOAMorg::triSurface
    OpenFOAMorg::twoPhaseModels_compressibleCavitation
    OpenFOAMorg::twoPhaseModels_compressibleInterfaceProperties
    OpenFOAMorg::twoPhaseModels_compressibleTwoPhases
    OpenFOAMorg::twoPhaseModels_incompressibleCavitation
    OpenFOAMorg::twoPhaseModels_incompressibleTwoPhases
    OpenFOAMorg::twoPhaseModels_interfaceCompression
    OpenFOAMorg::twoPhaseModels_interfaceProperties
    OpenFOAMorg::twoPhaseModels_twoPhaseMixture
    OpenFOAMorg::twoPhaseModels_twoPhaseProperties
    OpenFOAMorg::waves
)
# module vars
set(OpenFOAMorg_FOUND TRUE)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenFOAMorg
  FOUND_VAR OpenFOAMorg_FOUND
  REQUIRED_VARS OpenFOAMorg_FOUND
  VERSION_VAR $ENV{FOAM_API}
  HANDLE_COMPONENTS
)
