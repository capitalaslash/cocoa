import os
from pathlib import Path
import sys

from component import Component


COMPONENT_EXCLUDE = [
    "OpenFOAM",  # uses $(OBJECT_DIR)
    "OSspecific_POSIX",  # no library
    "Pstream_dummy",  # "OS specific"
    "Pstream_mpi",  # "OS specific"
    # "dummyThirdParty_MGridGen",  # external lib
    # "dummyThirdParty_metisDecomp",  # external lib
    # "dummyThirdParty_ptscotchDecomp",  # external lib
    # "dummyThirdParty_scotchDecomp",  # external lib
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

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} [org|com]")
        sys.exit(1)
    else:
        ver = sys.argv[1]
    print(f"generating cmake file for OpenFOAM.{ver}")
    Component.ver = ver

    try:
        of_dir = Path(os.environ["WM_PROJECT_DIR"])
        print(f"of_dir: {of_dir}")
    except KeyError:
        print("OpenFOAM directory not available")
        exit(1)

    # set libbin path
    path_libbin = os.environ["FOAM_LIBBIN"].replace(
        os.environ["WM_PROJECT_DIR"], f"${{OpenFOAM{ver}_DIR}}"
    )
    Component.path_libbin = path_libbin

    components = []
    lib_registry = {}

    # custom OpenFOAM base library
    component_of = Component(
        name="OpenFOAM",
        path=os.environ["WM_PROJECT_DIR"] + "/src/OpenFOAM",
        lib=os.environ["FOAM_LIBBIN"] + "/libOpenFOAM.so",
        includes=["${OpenFOAMorg_DIR}/src/OSspecific/POSIX/lnInclude"],
        dependencies=[],
        explicit_dependencies=[
            os.environ["FOAM_LIBBIN"] + "/openmpi-system/libPstream.so"
        ],
    )
    components.append(component_of)
    lib_registry["OpenFOAM"] = "OpenFOAM"

    # component_metis = Component(
    #     name="parallel_decompose_metisDecomp",
    #     path=os.environ["WM_PROJECT_DIR"] + "/src/dummyThirdParty/metisDecomp",
    #     lib=os.environ["FOAM_LIBBIN"] + "/dummy/libmetisDecomp.so",
    #     includes=["${OpenFOAMorg_DIR}/src/parallel/decompose/metisDecomp/lnInclude"],
    #     dependencies=[],
    # )
    # components.append(component_metis)
    # lib_registry["parallel_decompose_metisDecomp"] = "parallel_decompose_metisDecomp"

    for subdir, dirs, files in os.walk(of_dir / "src"):
        if (Path(subdir) / "Make").exists():
            trim_length = len(str(of_dir)) + 5
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
    with open(f"FindOpenFOAM{ver}.cmake", "w") as f:
        f.write(
            f"""
# OpenFOAM{ver} CMake module

if(NOT DEFINED OpenFOAM{ver}_DIR)
  if(EXISTS $ENV{{WM_PROJECT_DIR}})
    message(STATUS "using ENV{{WM_PROJECT_DIR}}")
    set(OpenFOAM{ver}_DIR "$ENV{{WM_PROJECT_DIR}}")
  else()
    message(
      FATAL_ERROR
      "OpenFOAM{ver}: no path available, set OpenFOAM{ver}_DIR or ENV{{WM_PROJECT_DIR}}"
    )
  endif()
endif()
message(STATUS "OpenFOAM{ver}_DIR: ${{OpenFOAM{ver}_DIR}}")

"""
        )

        f.write("# define components\n")
        [f.write(c.write_cmake()) for c in components]

        f.write("# set up dependencies\n")
        [f.write(c.write_deps_cmake()) for c in components]

        f.write(
            f"""
# definitions
target_compile_definitions(OpenFOAM{ver}::OpenFOAM
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
        f.write(f"add_library(OpenFOAM{ver}::ALL INTERFACE IMPORTED)\n")
        f.write(f"target_link_libraries(OpenFOAM{ver}::ALL\n  INTERFACE\n")
        [f.write(f"    OpenFOAM{ver}::{c.name}\n") for c in components]
        f.write(")\n")

        f.write(
            f"""
# module vars
set(OpenFOAM{ver}_FOUND TRUE)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenFOAM{ver}
  FOUND_VAR OpenFOAM{ver}_FOUND
  REQUIRED_VARS OpenFOAM{ver}_FOUND
  VERSION_VAR {os.environ["WM_PROJECT_VERSION"]}
  HANDLE_COMPONENTS
)
"""
        )
