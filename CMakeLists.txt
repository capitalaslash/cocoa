cmake_minimum_required(VERSION 3.15...3.29)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "type of build")
endif()

project(cocoa
  VERSION 0.1
  LANGUAGES CXX
)

message(STATUS "Version: ${PROJECT_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

set(LIB_FILES
  src/coupling_factory.hpp
  src/coupling_manager.hpp
  src/coupling_simple.hpp
  src/coupling_simple.cpp
  src/field_factory.hpp
  src/field_coupling.hpp
  src/problem_factory.hpp
  src/problem.hpp
  src/problem.cpp
  src/problem_fd1d.hpp
  src/problem_fd1d.cpp
  src/problem_femus.hpp
  src/problem_oforg.hpp
)
set(LIB_DEPENDENCIES
  fmt::fmt
)

# external dependencies ================================================

include(FeatureSummary)
include(FetchContent)

option(COCOA_FETCH_FMT "Get libfmt from the web" OFF)
if (COCOA_FETCH_FMT)
  FetchContent_Declare(fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG master
  )
  FetchContent_MakeAvailable(fmt)
else()
  find_package(fmt REQUIRED)
endif()

option(COCOA_ENABLE_OFORG "Enable OpenFOAM (org version)." OFF)
if (COCOA_ENABLE_OFCOM)
  find_package(OpenFOAMorg REQUIRED)
  list(APPEND LIB_FILES
    src/problem_oforg.hpp
    # src/problem_oforg.cpp
  )
  list(APPEND LIB_DEPENDENCIES
    OpenFOAMorg::ALL
  )
endif()

option(COCOA_ENABLE_PROXPDE "Enable ProXPDE library." OFF)
if (COCOA_ENABLE_PROXPDE)
  find_package(ProXPDE REQUIRED)
  list(APPEND LIB_FILES
    src/problem_proxpde.hpp
    src/problem_proxpde.cpp
  )
  list(APPEND LIB_DEPENDENCIES
    ProXPDE::proxpde
  )
endif()

option(COCOA_ENABLE_MEDCOUPLING "Enable MEDCoupling coupling." OFF)
if (COCOA_ENABLE_MEDCOUPLING)
  find_package(MEDCoupling
    REQUIRED
    PATHS
      $ENV{SPACK_ENV}/.spack-env/view/cmake_files
  )
  add_library(medcoupling INTERFACE)
  # fix missing target in medcoupling cmake
  target_link_libraries(medcoupling
    INTERFACE
      interpkernel
      medcouplingcpp
      medcouplingremapper
      # medicoco
      medloader
  )
  target_include_directories(medcoupling
    INTERFACE
      ${MEDCOUPLING_INCLUDE_DIRS}
  )
  add_library(medcoupling::medcoupling ALIAS medcoupling)
  list(APPEND LIB_FILES
    src/coupling_med.hpp
    src/coupling_med.cpp
    src/field_med.hpp
    src/field_med.cpp
    src/mesh_med.hpp
    src/mesh_med.cpp
  )
  list(APPEND LIB_DEPENDENCIES
    medcoupling::medcoupling
  )
endif()

find_package(Python 3.11 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# library ==============================================================

add_library(cocoa SHARED)
target_sources(cocoa
  PRIVATE
    ${LIB_FILES}
)
target_link_libraries(cocoa
  PUBLIC
  ${LIB_DEPENDENCIES}
)
add_library(cocoa::cocoa ALIAS cocoa)

configure_file(src/plugins.hpp.in plugins.hpp @ONLY)
target_include_directories(cocoa PUBLIC "${CMAKE_CURRENT_BINARY_DIR}")

# executable ===========================================================

add_executable(main)

set(MAIN_FILES
  src/main.cpp
)

target_sources(main
  PRIVATE
    ${MAIN_FILES}
	  # ${HDR_FILES}
)
target_link_libraries(main
  PUBLIC
    cocoa::cocoa
)

feature_summary(WHAT ALL)

# python ===============================================================

option(COCOA_ENABLE_PYTHON "Enable python bindings." OFF)
if (COCOA_ENABLE_PYTHON)
  pybind11_add_module(fd1d
    src/py_problem_fd1d.cpp
  )
  target_link_libraries(fd1d
    PRIVATE
      cocoa::cocoa
  )
  install(TARGETS fd1d DESTINATION .)

  pybind11_add_module(coupling_simple
    src/py_coupling_simple.cpp
  )
  target_link_libraries(coupling_simple  PRIVATE
      cocoa::cocoa
  )
  install(TARGETS coupling_simple DESTINATION .)

  set(PY_FILES
    # src/py_problem.cpp
    src/py_problem_fd1d.cpp
    src/py_coupling_simple.cpp
  )
endif()

# format ===============================================================

find_program(CLANGFORMAT clang-format)
if (CLANGFORMAT)
  set(SRC_FILES
    ${LIB_FILES}
    ${MAIN_FILES}
    ${PY_FILES}
  )
  add_custom_target(format
    COMMAND ${CLANGFORMAT}
    --verbose -i
    ${SRC_FILES} ${HEADER_FILES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
else()
  add_custom_target(format
    COMMAND echo "clang-format could not be found"
  )
endif()
