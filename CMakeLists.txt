cmake_minimum_required(VERSION 3.21)
project(GeoConstraint VERSION 1.0.0 LANGUAGES CXX
        DESCRIPTION "3D geometry and constraint-solving engine")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(GEOCONSTRAINT_BUILD_TESTS "Build the Catch2 test suite" ON)
option(GEOCONSTRAINT_BUILD_EXAMPLES "Build the example programs" ON)
option(GEOCONSTRAINT_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

# ---------------------------------------------------------------------------
# Dependencies. A system install is preferred; otherwise a pinned release
# archive (with SHA-256 verification) is downloaded at configure time.
# ---------------------------------------------------------------------------
include(FetchContent)

find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  message(STATUS "Eigen3 not found on the system; fetching pinned Eigen 3.4.0")
  FetchContent_Declare(
    eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_SUBDIR geoconstraint-no-cmake  # populate only; Eigen is header-only
  )
  FetchContent_MakeAvailable(eigen)
  add_library(geoconstraint_eigen INTERFACE)
  target_include_directories(geoconstraint_eigen SYSTEM INTERFACE "${eigen_SOURCE_DIR}")
  add_library(Eigen3::Eigen ALIAS geoconstraint_eigen)
endif()

# ---------------------------------------------------------------------------
# Warnings
# ---------------------------------------------------------------------------
add_library(geoconstraint_warnings INTERFACE)
if(MSVC)
  target_compile_options(geoconstraint_warnings INTERFACE /W4 /permissive-)
  if(GEOCONSTRAINT_WARNINGS_AS_ERRORS)
    target_compile_options(geoconstraint_warnings INTERFACE /WX)
  endif()
else()
  target_compile_options(geoconstraint_warnings INTERFACE
    -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
    -Wnon-virtual-dtor -Woverloaded-virtual -Wdouble-promotion)
  if(GEOCONSTRAINT_WARNINGS_AS_ERRORS)
    target_compile_options(geoconstraint_warnings INTERFACE -Werror)
  endif()
endif()

# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------
add_library(geoconstraint STATIC
  src/geometry/line_plane.cpp
  src/geometry/distances.cpp
  src/constraints/constraints.cpp
  src/solver/problem.cpp
  src/solver/diagnostics.cpp
  src/solver/solver.cpp)
add_library(GeoConstraint::geoconstraint ALIAS geoconstraint)
target_include_directories(geoconstraint PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(geoconstraint PUBLIC Eigen3::Eigen
                      PRIVATE geoconstraint_warnings)
target_compile_features(geoconstraint PUBLIC cxx_std_20)

enable_testing()

if(GEOCONSTRAINT_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
if(GEOCONSTRAINT_BUILD_TESTS)
  add_subdirectory(tests)
endif()
