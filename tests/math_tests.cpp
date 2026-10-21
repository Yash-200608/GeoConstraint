#include <cmath>
#include <numbers>

#include "geoconstraint/math/matrix3.hpp"
#include "geoconstraint/math/point3.hpp"
#include "geoconstraint/math/rigid_transform3.hpp"
#include "geoconstraint/math/vector3.hpp"
#include "test_utils.hpp"

using namespace geoconstraint;
using Catch::Matchers::WithinAbs;

namespace {
constexpr double kPi = std::numbers::pi;
}

TEST_CASE("Vector3 arithmetic with known values", "[math][vector]") {
  const Vector3 a{1, 2, 3};
  const Vector3 b{4, 5, 6};
  REQUIRE(a + b == Vector3{5, 7, 9});
  REQUIRE(b - a == Vector3{3, 3, 3});
  REQUIRE(-a == Vector3{-1, -2, -3});
  REQUIRE(a * 2.0 == Vector3{2, 4, 6});
  REQUIRE(2.0 * a == Vector3{2, 4, 6});
  REQUIRE(b / 2.0 == Vector3{2, 2.5, 3});
  REQUIRE(dot(a, b) == 32.0);
  REQUIRE(cross(a, b) == Vector3{-3, 6, -3});
  REQUIRE(squaredNorm(Vector3{3, 4, 12}) == 169.0);
  REQUIRE(norm(Vector3{3, 4, 12}) == 13.0);
}

TEST_CASE("Vector3 division by zero or non-finite scalar throws", "[math][vector][errors]") {
  const Vector3 v{1, 1, 1};
  REQUIRE_THROWS_AS(v / 0.0, InvalidArgumentError);
  REQUIRE_THROWS_AS(v / std::numeric_limits<double>::infinity(), InvalidArgumentError);
  REQUIRE_THROWS_AS(v / std::numeric_limits<double>::quiet_NaN(), InvalidArgumentError);
}

TEST_CASE("Cross product identities", "[math][vector]") {
  const Vector3 samples[] = {{1, 2, 3},     {-4, 0.5, 2},  {0, 0, 7},     {1e-3, -2e-3, 5e-4},
                             {100, -50, 25}, {0.3, 0.3, -0.9}, {-7, 8, -9}, {2, 2, 2}};
  for (const Vector3& a : samples) {
    for (const Vector3& b : samples) {
      const Vector3 c = cross(a, b);
      const double scale = norm(a) * norm(b);
      INFO("a = " << test::str(a) << " b = " << test::str(b));
      // c is orthogonal to both inputs
      REQUIRE(std::abs(dot(c, a)) <= 1e-12 * scale * norm(a));
      REQUIRE(std::abs(dot(c, b)) <= 1e-12 * scale * norm(b));
      // anticommutativity
      REQUIRE(approxEqual(c, -cross(b, a)));
      // Lagrange identity: |a x b|^2 + (a.b)^2 = |a|^2 |b|^2
      REQUIRE_NEAR(squaredNorm(c) + dot(a, b) * dot(a, b), squaredNorm(a) * squaredNorm(b),
                   1e-12 * squaredNorm(a) * squaredNorm(b) + 1e-300);
    }
  }
  // right-handedness
  REQUIRE(cross({1, 0, 0}, {0, 1, 0}) == Vector3{0, 0, 1});
  REQUIRE(cross({0, 1, 0}, {0, 0, 1}) == Vector3{1, 0, 0});
}

TEST_CASE("Normalization is safe", "[math][vector][errors]") {
  const Vector3 n = normalized({0, 3, 4});
  REQUIRE_NEAR(n.y, 0.6, 1e-15);
  REQUIRE_NEAR(n.z, 0.8, 1e-15);
  REQUIRE_NEAR(norm(n), 1.0, 1e-15);

  REQUIRE_THROWS_AS(normalized({0, 0, 0}), DegenerateGeometryError);
  REQUIRE_THROWS_AS(normalized({0, 1e-13, 0}), DegenerateGeometryError);  // below zeroLength
  REQUIRE_THROWS_AS(normalized({std::nan(""), 0, 1}), DegenerateGeometryError);
  REQUIRE_THROWS_AS(normalized({std::numeric_limits<double>::infinity(), 0, 0}),
                    DegenerateGeometryError);
  REQUIRE_FALSE(tryNormalized({0, 0, 0}).has_value());
  REQUIRE(tryNormalized({0, 0, 2}).has_value());
}

TEST_CASE("Angles between vectors", "[math][vector]") {
  REQUIRE_NEAR(angleBetween({1, 0, 0}, {0, 1, 0}), kPi / 2, 1e-15);
  REQUIRE_NEAR(angleBetween({1, 0, 0}, {1, 0, 0}), 0.0, 1e-15);
  REQUIRE_NEAR(angleBetween({1, 0, 0}, {-3, 0, 0}), kPi, 1e-15);
  REQUIRE_NEAR(angleBetween({2, 0, 0}, {1, std::sqrt(3.0), 0}), kPi / 3, 1e-15);
  REQUIRE_NEAR(angleBetween({1, 1, 0}, {0, 1, 0}), kPi / 4, 1e-15);
  REQUIRE_THROWS_AS(angleBetween({0, 0, 0}, {1, 0, 0}), DegenerateGeometryError);
  REQUIRE_THROWS_AS(angleBetween({1, 0, 0}, {0, 0, 0}), DegenerateGeometryError);
}

TEST_CASE("Parallel and perpendicular predicates", "[math][vector]") {
  REQUIRE(isParallel({1, 2, 3}, {2, 4, 6}));
  REQUIRE(isParallel({1, 2, 3}, {-2, -4, -6}));  // antiparallel counts as parallel
  REQUIRE_FALSE(isParallel({1, 0, 0}, {0, 1, 0}));
  REQUIRE(isPerpendicular({1, 0, 0}, {0, 5, 0}));
  REQUIRE(isPerpendicular({1, 1, 0}, {1, -1, 5}));
  REQUIRE_FALSE(isPerpendicular({1, 0, 0}, {1, 1, 0}));
  // a zero vector is never "parallel" or "perpendicular" - it is an error
  REQUIRE_THROWS_AS(isParallel({0, 0, 0}, {1, 0, 0}), DegenerateGeometryError);
  REQUIRE_THROWS_AS(isPerpendicular({0, 0, 0}, {1, 0, 0}), DegenerateGeometryError);
}

TEST_CASE("Approximate equality", "[math][vector][tolerances]") {
  REQUIRE(approxEqual(Vector3{1, 2, 3}, Vector3{1, 2, 3 + 1e-10}));
  REQUIRE_FALSE(approxEqual(Vector3{1, 2, 3}, Vector3{1, 2, 3 + 1e-6}));
  REQUIRE_FALSE(approxEqual(Vector3{1, 2, std::nan("")}, Vector3{1, 2, std::nan("")}));
  REQUIRE(approxEqual(1.0, 1.0 + 5e-10));
  REQUIRE_FALSE(approxEqual(1.0, 1.0 + 1e-7));
  // tolerances can be tightened explicitly
  Tolerances strict;
  strict.length = 1e-14;
  strict.relative = 1e-14;
  REQUIRE_FALSE(approxEqual(Vector3{1, 2, 3}, Vector3{1, 2, 3 + 1e-10}, strict));
}

TEST_CASE("Point3 and Vector3 are distinct types", "[math][point]") {
  const Point3 p{1, 2, 3};
  const Point3 q{4, 6, 3};
  const Vector3 d = q - p;
  REQUIRE(d == Vector3{3, 4, 0});
  REQUIRE(p + d == q);
  REQUIRE(q - d == p);
  REQUIRE(distance(p, q) == 5.0);
  REQUIRE(squaredDistance(p, q) == 25.0);
  REQUIRE(approxEqual(p, Point3{1, 2, 3 + 1e-12}));
  REQUIRE_FALSE(approxEqual(p, q));
  static_assert(!std::is_convertible_v<Point3, Vector3>);
  static_assert(!std::is_convertible_v<Vector3, Point3>);
}

TEST_CASE("Matrix3 basics and rotations", "[math][matrix]") {
  const Matrix3 rz = Matrix3::rotationAxisAngle({0, 0, 1}, kPi / 2);
  const Vector3 y = rz * Vector3{1, 0, 0};
  REQUIRE_NEAR(y.x, 0.0, 1e-15);
  REQUIRE_NEAR(y.y, 1.0, 1e-15);
  REQUIRE_NEAR(y.z, 0.0, 1e-15);
  REQUIRE(rz.isRotation());
  REQUIRE_NEAR(rz.determinant(), 1.0, 1e-15);
  REQUIRE(Matrix3::identity().isRotation());

  // arbitrary axis: R^T R = I and the axis is fixed
  const Vector3 axis{1, 2, -2};
  const Matrix3 r = Matrix3::rotationAxisAngle(axis, 1.234);
  REQUIRE(r.isRotation());
  REQUIRE(approxEqual(r * axis, axis));
  const Matrix3 id = r.transposed() * r;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) REQUIRE_NEAR(id(i, j), i == j ? 1.0 : 0.0, 1e-14);
  }

  // Rotation by angle preserves lengths and the angle between vectors.
  const Vector3 a{1, -2, 0.5};
  const Vector3 b{0.3, 4, -1};
  REQUIRE_NEAR(norm(r * a), norm(a), 1e-14);
  REQUIRE_NEAR(angleBetween(r * a, r * b), angleBetween(a, b), 1e-14);
}

TEST_CASE("Matrix3 rejects non-rotations and bad rotation arguments", "[math][matrix][errors]") {
  const Matrix3 reflection{{-1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  const Matrix3 scaling{{2, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  REQUIRE_FALSE(reflection.isRotation());  // det = -1
  REQUIRE_FALSE(scaling.isRotation());     // R^T R != I
  REQUIRE_THROWS_AS(Matrix3::rotationAxisAngle({0, 0, 0}, 1.0), DegenerateGeometryError);
  REQUIRE_THROWS_AS(Matrix3::rotationAxisAngle({0, 0, 1}, std::nan("")), InvalidArgumentError);
}

TEST_CASE("RigidTransform3 applies to points and vectors differently", "[math][transform]") {
  const RigidTransform3 t = RigidTransform3::translation({1, 2, 3});
  REQUIRE(t.apply(Point3{1, 1, 1}) == Point3{2, 3, 4});
  REQUIRE(t.apply(Vector3{1, 1, 1}) == Vector3{1, 1, 1});  // directions ignore translation

  const RigidTransform3 r = RigidTransform3::rotation({0, 0, 1}, kPi / 2);
  const Point3 rotated = r.apply(Point3{1, 0, 0});
  REQUIRE_POINT_NEAR(rotated, (Point3{0, 1, 0}), 1e-15);
}

TEST_CASE("RigidTransform3 composition order and inverse", "[math][transform]") {
  const RigidTransform3 rot = RigidTransform3::rotation({0, 0, 1}, kPi / 2);
  const RigidTransform3 shift = RigidTransform3::translation({10, 0, 0});

  // (shift * rot)(p) = shift(rot(p)): rotate first, then translate
  const Point3 a = (shift * rot).apply(Point3{1, 0, 0});
  REQUIRE_POINT_NEAR(a, (Point3{10, 1, 0}), 1e-14);
  // (rot * shift)(p) = rot(shift(p)): translate first, then rotate
  const Point3 b = (rot * shift).apply(Point3{1, 0, 0});
  REQUIRE_POINT_NEAR(b, (Point3{0, 11, 0}), 1e-14);

  // composition equals sequential application for a general chain
  const RigidTransform3 t1 = RigidTransform3::rotationAbout({1, 1, 0}, {1, 1, 1}, 0.7);
  const RigidTransform3 t2 = RigidTransform3::rotation({0, 1, 0}, -1.1);
  const RigidTransform3 t3 = RigidTransform3::translation({-2, 5, 0.5});
  const Point3 p{3, -1, 2};
  REQUIRE_POINT_NEAR(((t3 * t2) * t1).apply(p), t3.apply(t2.apply(t1.apply(p))), 1e-13);
  REQUIRE_POINT_NEAR((t3 * (t2 * t1)).apply(p), ((t3 * t2) * t1).apply(p), 1e-13);

  // inverse
  const RigidTransform3 all = (t3 * t2) * t1;
  REQUIRE_POINT_NEAR(all.inverse().apply(all.apply(p)), p, 1e-13);
  REQUIRE_POINT_NEAR((all * all.inverse()).apply(p), p, 1e-13);
}

TEST_CASE("RigidTransform3 preserves distances and angles", "[math][transform]") {
  const RigidTransform3 t = RigidTransform3::rotationAbout({4, -2, 1}, {2, -1, 3}, 2.2) *
                            RigidTransform3::translation({5, 5, 5});
  const Point3 p{1, 2, 3};
  const Point3 q{-4, 0.5, 7};
  const Point3 r{0, 9, -3};
  REQUIRE_NEAR(distance(t.apply(p), t.apply(q)), distance(p, q), 1e-13);
  REQUIRE_NEAR(angleBetween(t.apply(Vector3(q - p)), t.apply(Vector3(r - p))),
               angleBetween(q - p, r - p), 1e-13);
}

TEST_CASE("rotationAbout leaves the center fixed", "[math][transform]") {
  const Point3 center{3, 4, 5};
  const RigidTransform3 t = RigidTransform3::rotationAbout(center, {0, 1, 0}, kPi);
  REQUIRE_POINT_NEAR(t.apply(center), center, 1e-14);
  // rotating (4,4,5) by pi about the vertical axis through (3,4,5) gives (2,4,5)
  REQUIRE_POINT_NEAR(t.apply(Point3{4, 4, 5}), (Point3{2, 4, 5}), 1e-14);
}

TEST_CASE("RigidTransform3 rejects invalid input", "[math][transform][errors]") {
  const Matrix3 reflection{{-1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  REQUIRE_THROWS_AS(RigidTransform3(reflection, Vector3{}), InvalidArgumentError);
  REQUIRE_THROWS_AS(RigidTransform3(Matrix3::identity(), Vector3{std::nan(""), 0, 0}),
                    InvalidArgumentError);
}
