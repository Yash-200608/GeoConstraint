#pragma once

#include <cmath>
#include <sstream>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "geoconstraint/math/point3.hpp"

namespace test {

inline std::string str(const geoconstraint::Point3& p) {
  std::ostringstream os;
  os.precision(17);
  os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
  return os.str();
}

inline std::string str(const geoconstraint::Vector3& v) {
  return str(geoconstraint::Point3{v.x, v.y, v.z});
}

}  // namespace test

#define REQUIRE_NEAR(actual, expected, tol) \
  REQUIRE_THAT((actual), Catch::Matchers::WithinAbs((expected), (tol)))

// Euclidean distance between two points/vectors must be <= tol (absolute).
#define REQUIRE_POINT_NEAR(actual, expected, tol)                                  \
  do {                                                                             \
    INFO("actual = " << test::str(actual) << ", expected = " << test::str(expected)); \
    REQUIRE(geoconstraint::norm(geoconstraint::Vector3{(actual).x - (expected).x,  \
                                                       (actual).y - (expected).y,  \
                                                       (actual).z - (expected).z}) <= (tol)); \
  } while (false)
