import os
from pathlib import Path

from component import Component


COMPONENT_EXCLUDE = [
    "OpenFOAM",  # uses $(OBJECT_DIR)
    "OSspecific_POSIX",  # no library
    "Pstream_dummy",  # "OS specific"
    "Pstream_mpi",  # "OS specific"
    "dummyThirdParty_MGridGen",  # external lib
    "dummyThirdParty_metisDecomp",  # external lib
    "dummyThirdParty_ptscotchDecomp",  # external lib
    "dummyThirdParty_scotchDecomp",  # external lib
    "fvAgglomerationMethods_MGridGenGamgAgglomeration",  # disabled by default
    "fvMeshDistributors",  # depends on optional libs
    "parallel_decompose_metisDecomp",  # external lib
    "parallel_decompose_ptscotchDecomp",  # external lib
    "parallel_decompose_scotchDecomp",  # external lib
    "parallel_decompose_zoltanDecomp",  # external lib
    "renumber_SloanRenumber",  # external lib
    "renumber_zoltanRenumber",  # external lib
    # "thermophysicalModels_chemistryModel",  # missing include
]

try:
    oforg_dir = Path(os.environ["WM_PROJECT_DIR"])
    print(f"oforg_dir: {oforg_dir}")
except KeyError:
    print("OpenFOAM directory not available")
    exit(1)

# set libbin path
path_libbin = os.environ["FOAM_LIBBIN"].replace(
    os.environ["WM_PROJECT_DIR"], "${OpenFOAMorg_DIR}"
)
Component.path_libbin = path_libbin

components = []
lib_registry = {}

# OpenFOAM base library
base_component = Component(
    name="OpenFOAM",
    path=os.environ["WM_PROJECT_DIR"] + "/src/OpenFOAM",
    lib=os.environ["FOAM_LIBBIN"] + "/libOpenFOAM.so",
    includes=["${OpenFOAMorg_DIR}/src/OSspecific/POSIX/lnInclude"],
    dependencies=[],
    explicit_dependencies=[os.environ["FOAM_LIBBIN"] + "/openmpi-system/libPstream.so"],
)
components.append(base_component)
lib_registry["OpenFOAM"] = "OpenFOAM"

for subdir, dirs, files in os.walk(oforg_dir / "src"):
    if (Path(subdir) / "Make").exists():
        trim_length = len(str(oforg_dir)) + 5
        name = subdir[trim_length:].replace("/", "_")
        if name not in COMPONENT_EXCLUDE:
            new_component = Component.from_path(name, subdir)
            components.append(new_component)
            libname = os.path.basename(new_component.lib)[3:-3]
            lib_registry[libname] = new_component.name
        else:
            print(f"{name} excluded")
[c.sanitize_dependencies(lib_registry) for c in components]

# print out the cmake module
with open("FindOpenFOAMorg.cmake", "w") as f:
    f.write(
        """
# OpenFOAMorg CMake module

if(NOT DEFINED OpenFOAMorg_DIR)
  if(EXISTS $ENV{WM_PROJECT_DIR})
    message(STATUS "using ENV{WM_PROJECT_DIR}")
    set(OpenFOAMorg_DIR "$ENV{WM_PROJECT_DIR}")
  else()
    message(FATAL_ERROR "OpenFOAMorg: no path available, set OpenFOAMorg_DIR or ENV{WM_PROJECT_DIR}")
  endif()
endif()
message(STATUS "OpenFOAMorg_DIR: ${OpenFOAMorg_DIR}")

"""
    )

    f.write("# define components\n")
    [f.write(c.write_cmake()) for c in components]

    f.write("# set up dependencies\n")
    [f.write(c.write_deps_cmake()) for c in components]

    f.write(
        """
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

"""
    )

    f.write("# add ::ALL target\n")
    f.write("add_library(OpenFOAMorg::ALL INTERFACE IMPORTED)\n")
    f.write("target_link_libraries(OpenFOAMorg::ALL\n  INTERFACE\n")
    [f.write(f"    OpenFOAMorg::{c.name}\n") for c in components]
    f.write(")\n")

    f.write(
        f"""
# module vars
set(OpenFOAMorg_FOUND TRUE)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenFOAMorg
  FOUND_VAR OpenFOAMorg_FOUND
  REQUIRED_VARS OpenFOAMorg_FOUND
  VERSION_VAR {os.environ["WM_PROJECT_VERSION"]}
  HANDLE_COMPONENTS
)
"""
    )
