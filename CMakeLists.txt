cmake_minimum_required(VERSION 3.25...3.30)

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

set(LIB_SOURCES
  src/coupling/coupling_simple.cpp
  src/coupling/coupling_factory.cpp
  src/coupling/field_factory.cpp
  src/coupling/mesh_factory.cpp
  src/problem/fdutils.cpp
  src/problem/problem_fd1d.cpp
  src/problem/problem_fd2d.cpp
  src/problem/problem_factory.cpp
)
set(LIB_HEADERS
  src/enums.hpp
  src/coupling/coupling_manager.hpp
  src/coupling/coupling_simple.hpp
  src/coupling/field_coupling.hpp
  src/coupling/mesh_coupling.hpp
  src/problem/problem.hpp
  src/problem/fdutils.hpp
  src/problem/problem_fd1d.hpp
  src/problem/problem_fd2d.hpp
  src/problem/problem_femus.hpp
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
  list(APPEND LIB_SOURCES
    src/problem/problem_oforg.cpp
  )
  list(APPEND LIB_HEADERS
    src/problem/problem_oforg.hpp
  )
  list(APPEND LIB_DEPENDENCIES
    OpenFOAMorg::ALL
    $ENV{FOAM_LIBBIN}/dummy/libmetisDecomp.so
  )
endif()

# proxpde ---
if (COCOA_ENABLE_PROXPDE)
  find_package(ProXPDE REQUIRED)
  list(APPEND LIB_SOURCES
    src/problem/problem_proxpde.cpp
  )
  list(APPEND LIB_HEADERS
    src/problem/problem_proxpde.hpp
  )
  list(APPEND LIB_DEPENDENCIES
    ProXPDE::proxpde
  )
endif()

# medcoupling ---
if (COCOA_ENABLE_MEDCOUPLING)
  find_package(MEDCoupling REQUIRED
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
  # target_include_directories(medcoupling SYSTEM
  #   INTERFACE
  #     ${MEDCOUPLING_INCLUDE_DIRS}
  # )
  add_library(medcoupling::medcoupling ALIAS medcoupling)
  list(APPEND LIB_SOURCES
    src/coupling/coupling_med.cpp
    src/coupling/field_med.cpp
    src/coupling/mesh_med.cpp
  )
  list(APPEND LIB_HEADERS
    src/coupling/coupling_med.hpp
    src/coupling/field_med.hpp
    src/coupling/mesh_med.hpp
  )
  list(APPEND LIB_DEPENDENCIES
    medcoupling::medcoupling
  )
endif()

# pybind11 ---
if (COCOA_ENABLE_PYTHON)
  find_package(Python 3.10 COMPONENTS Interpreter Development REQUIRED)
  find_package(pybind11 CONFIG REQUIRED)
  find_package(Eigen3 3.4 REQUIRED)
endif()


# library ==============================================================

add_library(cocoa)
target_sources(cocoa
  PRIVATE
    ${LIB_SOURCES}
    ${CMAKE_CURRENT_BINARY_DIR}/plugins.hpp
)
target_sources(cocoa PUBLIC
  FILE_SET HEADERS
  BASE_DIRS src
  FILES ${LIB_HEADERS}
)
target_link_libraries(cocoa
  PUBLIC
    ${LIB_DEPENDENCIES}
)
target_compile_features(cocoa PRIVATE cxx_std_20)
add_library(cocoa::cocoa ALIAS cocoa)

configure_file(src/plugins.hpp.in plugins.hpp @ONLY)
target_include_directories(cocoa
  PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}
)

feature_summary(WHAT ALL)


# python ===============================================================

if (COCOA_ENABLE_PYTHON)
  pybind11_add_module(pycocoa
    python/bindings.cpp
  )
  target_link_libraries(pycocoa
    PRIVATE
      cocoa::cocoa
      Eigen3::Eigen
  )
  target_include_directories(pycocoa
    PRIVATE
      ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()


# test =================================================================

if (COCOA_ENABLE_TESTS)
  include(CTest)

  if (COCOA_ENABLE_PYTHON)
    macro(cocoa_add_test_python name)
      message(STATUS "adding python app ${name}")

      add_test(
        NAME cocoa::python::${name}
        COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/python/app/${name}run.py
      )
      set_property(TEST cocoa::python::${name}
        PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}
      )

      set(extra_args ${ARGN})
      list(LENGTH extra_args extra_count)
      if (${extra_count} GREATER 0)
        foreach(datafile ${extra_args})
          MESSAGE(STATUS "copying ${datafile} to binary dir")
          file(
            COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/${datafile}
            DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
          )
        endforeach()
      endif()

    endmacro()

    cocoa_add_test_python(fd1d fd1d_heat.dat)
    cocoa_add_test_python(couplingsimple fd1d_heat.dat fd1d_hc.dat)
    cocoa_add_test_python(fd2d fd2d_heat.dat)
    cocoa_add_test_python(oforg oforg1 oforg1.dat)
    cocoa_add_test_python(proxpde proxpde_heat.yaml)
    cocoa_add_test_python(couplingmed proxpde_heat.yaml proxpde_hc.yaml)
    cocoa_add_test_python(buoyant proxpde_buoyant_ns.yaml proxpde_buoyant_heat.yaml)
    cocoa_add_test_python(custom)

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
  target_include_directories(app_${name}
    PRIVATE
      ${CMAKE_CURRENT_BINARY_DIR}
  )

  if (COCOA_ENABLE_TESTS)
    add_test(NAME cocoa::app::${name} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/app_${name})
  endif()

  set(extra_args ${ARGN})
  list(LENGTH extra_args extra_count)
  if (${extra_count} GREATER 0)
    foreach(datafile ${extra_args})
      MESSAGE(STATUS "copying ${datafile} to binary dir")
      file(
        COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/${datafile}
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
      )
    endforeach()
  endif()

endmacro()

cocoa_add_app(fd1d fd1d_heat.dat)
cocoa_add_app(fd2d fd2d_heat.dat)
cocoa_add_app(simple fd1d_heat.dat fd1d_hc.dat)
cocoa_add_app(cyclic fd2d_heatpipe.dat)

if (COCOA_ENABLE_PROXPDE)
  cocoa_add_app(proxpde proxpde_heat.yaml proxpde_hc.yaml)
endif()

if (COCOA_ENABLE_OFORG)
  cocoa_add_app(oforg oforg1)
endif()

cocoa_add_app(ns proxpde_ns.yaml)
cocoa_add_app(buoyant proxpde_buoyant_ns.yaml proxpde_buoyant_heat.yaml)


# install ==============================================================

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cocoa-config-version.cmake
  COMPATIBILITY AnyNewerVersion
)

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/cocoa-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cocoa
)

install(
  TARGETS cocoa
  EXPORT cocoaTargets
  FILE_SET HEADERS
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cocoa
)

if (COCOA_ENABLE_MEDCOUPLING)
  install(
    TARGETS medcoupling
    EXPORT cocoaTargets
  )
endif()

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/plugins.hpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cocoa
)

if (COCOA_ENABLE_PYTHON)
  install(
    TARGETS pycocoa
    EXPORT cocoaTargets
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
endif()

install(
  EXPORT cocoaTargets
  NAMESPACE cocoa::
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cocoa
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/cocoa-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cocoa-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_DATADIR}/cocoa
)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/cocoa-config.cmake
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cocoa
)

install(
  FILES ${CMAKE_CURRENT_LIST_DIR}/cmake/FindOpenFOAMorg.cmake
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cocoa
)


# format ===============================================================

find_program(CLANGFORMAT clang-format)
if (CLANGFORMAT)
  set(SRC_FILES
    ${LIB_SOURCES}
    ${LIB_HEADERS}
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
