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
  src/coupling_manager.hpp
  src/coupling_simple.hpp
  src/coupling_simple.cpp
  src/coupling_factory.cpp
  src/enums.hpp
  src/field_coupling.hpp
  src/field_factory.cpp
  src/fdutils.hpp
  src/fdutils.cpp
  src/mesh_coupling.hpp
  src/mesh_factory.cpp
  src/problem.hpp
  src/problem_fd1d.hpp
  src/problem_fd1d.cpp
  src/problem_femus.hpp
  src/problem_oforg.hpp
  src/problem_factory.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/plugins.hpp
)
set(LIB_DEPENDENCIES
  fmt::fmt
)

# options ==============================================================

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)
option(COCOA_ENABLE_OFORG "Enable OpenFOAM (org version)." OFF)
option(COCOA_ENABLE_MEDCOUPLING "Enable MEDCoupling coupling." OFF)
option(COCOA_ENABLE_PROXPDE "Enable ProXPDE library." OFF)
option(COCOA_ENABLE_PYTHON "Enable python bindings." OFF)
option(COCOA_ENABLE_TESTS "Enable tests." ON)
option(COCOA_FETCH_FMT "Get libfmt from the web" OFF)

# external dependencies ================================================

include(FeatureSummary)
include(FetchContent)

# fmt ---
if (COCOA_FETCH_FMT)
  FetchContent_Declare(fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG master
  )
  FetchContent_MakeAvailable(fmt)
else()
  find_package(fmt REQUIRED)
endif()

# openfoam.org ---
if (COCOA_ENABLE_OFORG)
  find_package(OpenFOAMorg REQUIRED)
  list(APPEND LIB_FILES
    src/problem_oforg.hpp
    src/problem_oforg.cpp
  )
  list(APPEND LIB_DEPENDENCIES
    OpenFOAMorg::ALL
    $ENV{FOAM_LIBBIN}/dummy/libmetisDecomp.so
  )
endif()

# proxpde ---
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

# medcoupling ---
if (COCOA_ENABLE_MEDCOUPLING)
  find_package(MEDCoupling
    REQUIRED
    PATHS
      ${CMAKE_SOURCE_DIR}/contrib/medcoupling/cmake_files
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
  target_include_directories(medcoupling SYSTEM
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

# pybind11 ---
if (COCOA_ENABLE_PYTHON)
  find_package(Python 3.11 COMPONENTS Interpreter Development REQUIRED)
  find_package(pybind11 CONFIG REQUIRED)
endif()

# library ==============================================================

add_library(cocoa)
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
target_include_directories(cocoa
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}
)

feature_summary(WHAT ALL)

# python ===============================================================

if (COCOA_ENABLE_PYTHON)
  pybind11_add_module(pycocoa
    src/python_bindings.cpp
  )
  target_link_libraries(pycocoa
    PRIVATE
      cocoa::cocoa
  )
  install(TARGETS pycocoa DESTINATION .)
endif()

# test =================================================================

if (COCOA_ENABLE_TESTS)
  include(CTest)

  if (COCOA_ENABLE_PYTHON)
    macro(cocoa_add_test_python name)
      add_test(NAME cocoa::python::${name} COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${name}run.py)
      set_property(TEST cocoa::python::${name} PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR})
    endmacro()

    cocoa_add_test_python(fd1d)
    cocoa_add_test_python(couplingsimple)
    cocoa_add_test_python(proxpde)
    cocoa_add_test_python(couplingmed)
    cocoa_add_test_python(buoyant)

  endif()
endif()

# application ==========================================================

set(APP_FILES "")

macro(cocoa_add_app name)
  message(STATUS "adding app ${name}")
  add_executable(app_${name})
  list(APPEND APP_FILES app/app_${name}.cpp)
  target_sources(app_${name}
    PRIVATE
      app/app_${name}.cpp
  )
  target_link_libraries(app_${name}
    PUBLIC
      cocoa::cocoa
  )

  if (COCOA_ENABLE_TESTS)
    add_test(NAME cocoa::app::${name} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/app_${name})
  endif()
endmacro()

cocoa_add_app(fd1d)
cocoa_add_app(simple)

if (COCOA_ENABLE_PROXPDE)
  cocoa_add_app(proxpde)
endif()

if (COCOA_ENABLE_OFORG)
  cocoa_add_app(oforg)
endif()

cocoa_add_app(ns)
cocoa_add_app(buoyant)


# format ===============================================================

find_program(CLANGFORMAT clang-format)
if (CLANGFORMAT)
  set(SRC_FILES
    ${LIB_FILES}
    ${APP_FILES}
    src/python_bindings.cpp
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
