find_package(Catch2 3 QUIET)
if(NOT TARGET Catch2::Catch2WithMain)
  message(STATUS "Catch2 v3 not found on the system; fetching pinned Catch2 3.7.1")
  FetchContent_Declare(
    Catch2
    URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.7.1.tar.gz
    URL_HASH SHA256=c991b247a1a0d7bb9c39aa35faf0fe9e19764213f28ffba3109388e62ee0269c
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_MakeAvailable(Catch2)
  list(APPEND CMAKE_MODULE_PATH "${catch2_SOURCE_DIR}/extras")
endif()
include(Catch)

add_executable(geoconstraint_tests
  math_tests.cpp geometry_tests.cpp constraint_tests.cpp solver_tests.cpp numerical_tests.cpp)
target_link_libraries(geoconstraint_tests PRIVATE
  geoconstraint geoconstraint_warnings Catch2::Catch2WithMain)
catch_discover_tests(geoconstraint_tests DISCOVERY_MODE PRE_TEST)
