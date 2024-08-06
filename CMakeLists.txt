cmake_minimum_required(VERSION 3.15...3.29)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "type of build")
endif()

project(medhub
  VERSION 0.1
  LANGUAGES CXX
)

message(STATUS "Version: ${PROJECT_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic")

# list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# external dependencies ================================================

include(FeatureSummary)
include(FetchContent)

option(MEDHUB_FETCH_FMT "Get libfmt from the web" False)
if (MEDHUB_FETCH_FMT)
  FetchContent_Declare(fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG master
  )
  FetchContent_MakeAvailable(fmt)
else()
  find_package(fmt REQUIRED)
endif()

find_package(ProXPDE REQUIRED)

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

# library ==============================================================

add_library(medhub SHARED)
set(LIB_FILES
  src/med_manager.hpp
  src/med_field.hpp
  src/med_field.cpp
  src/med_mesh.hpp
  src/med_mesh.cpp
  src/problem.hpp
  src/problem_femus.hpp
  src/problem_oforg.hpp
  src/problem_proxpde.hpp
  src/problem_proxpde.cpp
  src/problem_fd1d.hpp
  src/problem_fd1d.cpp
)

target_sources(medhub
  PRIVATE
    ${LIB_FILES}
)
target_link_libraries(medhub
  PUBLIC
    medcoupling::medcoupling
    ProXPDE::proxpde
)
add_library(medhub::medhub ALIAS medhub)

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
    medhub::medhub
    fmt::fmt
)

feature_summary(WHAT ALL)

# format ===============================================================

find_program(CLANGFORMAT clang-format)
if (CLANGFORMAT)
  set(SRC_FILES
    ${LIB_FILES}
    ${MAIN_FILES}
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
