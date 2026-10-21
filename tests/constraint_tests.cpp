#include <cmath>
#include <memory>
#include <numbers>
#include <vector>

#include "geoconstraint/constraints/angle_constraint.hpp"
#include "geoconstraint/constraints/coincidence_constraint.hpp"
#include "geoconstraint/constraints/distance_constraint.hpp"
#include "geoconstraint/constraints/fixed_point_constraint.hpp"
#include "geoconstraint/constraints/parallel_constraint.hpp"
#include "geoconstraint/constraints/perpendicular_constraint.hpp"
#include "geoconstraint/constraints/point_on_plane_constraint.hpp"
#include "geoconstraint/solver/problem.hpp"
#include "geoconstraint/solver/solver.hpp"
#include "test_utils.hpp"

using namespace geoconstraint;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr PointId P0{0}, P1{1}, P2{2}, P3{3};

// Four points so that DirectionRef{P0,P1} and DirectionRef{P2,P3} are available.
std::vector<Point3> lines(const Vector3& u, const Vector3& v) {
  return {{0, 0, 0}, Point3{0, 0, 0} + u, {5, 5, 5}, Point3{5, 5, 5} + v};
}
const DirectionRef kDirA{P0, P1};
const DirectionRef kDirB{P2, P3};

}  // namespace

TEST_CASE("CoincidenceConstraint: r = p - q", "[constraint][coincident]") {
  const CoincidenceConstraint c(P0, P1);
  REQUIRE(c.typeName() == "Coincident");
  REQUIRE(c.residualCount() == 3);
  const std::vector<Point3> same = {{1, 2, 3}, {1, 2, 3}};
  const auto zero = c.residuals(same);
  REQUIRE(zero.size() == 3);
  for (double r : zero) REQUIRE(r == 0.0);

  const std::vector<Point3> apart = {{1, 2, 3}, {0, 2, 5}};
  const auto r = c.residuals(apart);
  REQUIRE(r == std::vector<double>{1, 0, -2});
}

TEST_CASE("FixedPointConstraint: r = p - p0", "[constraint][fixed]") {
  const FixedPointConstraint c(P0, Point3{1, 1, 1});
  REQUIRE(c.residualCount() == 3);
  REQUIRE(c.residuals(std::vector<Point3>{{1, 1, 1}}) == std::vector<double>{0, 0, 0});
  REQUIRE(c.residuals(std::vector<Point3>{{1, 2, 3}}) == std::vector<double>{0, 1, 2});
  REQUIRE_THROWS_AS(FixedPointConstraint(P0, Point3{std::nan(""), 0, 0}), InvalidArgumentError);
}

TEST_CASE("DistanceConstraint: r = |p - q| - d", "[constraint][distance]") {
  const std::vector<Point3> pts = {{0, 0, 0}, {3, 4, 0}};
  REQUIRE(DistanceConstraint(P0, P1, 5.0).residuals(pts)[0] == 0.0);   // satisfied
  REQUIRE(DistanceConstraint(P0, P1, 2.0).residuals(pts)[0] == 3.0);   // too far apart
  REQUIRE(DistanceConstraint(P0, P1, 7.0).residuals(pts)[0] == -2.0);  // too close
  REQUIRE(DistanceConstraint(P0, P1, 5.0).residualCount() == 1);
  // d = 0 is legal (equivalent to coincidence, but with a non-smooth residual)
  const std::vector<Point3> coincident = {{1, 1, 1}, {1, 1, 1}};
  REQUIRE(DistanceConstraint(P0, P1, 0.0).residuals(coincident)[0] == 0.0);

  REQUIRE_THROWS_AS(DistanceConstraint(P0, P1, -1.0), InvalidArgumentError);
  REQUIRE_THROWS_AS(DistanceConstraint(P0, P1, std::nan("")), InvalidArgumentError);
  REQUIRE_THROWS_AS(
      DistanceConstraint(P0, P1, std::numeric_limits<double>::infinity()), InvalidArgumentError);
}

TEST_CASE("PointOnPlaneConstraint: r = n . (p - p0)", "[constraint][plane]") {
  // plane z = 2; normal supplied non-unit and normalized by Plane3
  const PointOnPlaneConstraint c(P0, Plane3::fromPointNormal({0, 0, 2}, {0, 0, 3}));
  REQUIRE(c.residualCount() == 1);
  REQUIRE(c.residuals(std::vector<Point3>{{5, -5, 2}})[0] == 0.0);
  REQUIRE_NEAR(c.residuals(std::vector<Point3>{{0, 0, 5}})[0], 3.0, 1e-15);    // above
  REQUIRE_NEAR(c.residuals(std::vector<Point3>{{0, 0, -1}})[0], -3.0, 1e-15);  // below

  // plane x + y = 0 through the origin; (2,0,0) is sqrt(2) off the plane
  const PointOnPlaneConstraint tilted(P0, Plane3::fromPointNormal({0, 0, 0}, {1, 1, 0}));
  REQUIRE_NEAR(tilted.residuals(std::vector<Point3>{{2, 0, 0}})[0], std::sqrt(2.0), 1e-15);
  REQUIRE_NEAR(tilted.residuals(std::vector<Point3>{{1, -1, 9}})[0], 0.0, 1e-15);
}

TEST_CASE("ParallelConstraint: r = u x v", "[constraint][parallel]") {
  const ParallelConstraint c(kDirA, kDirB);
  REQUIRE(c.typeName() == "Parallel");
  REQUIRE(c.residualCount() == 3);

  // perpendicular unit directions: |u x v| = 1
  auto r = c.residuals(lines({1, 0, 0}, {0, 1, 0}));
  REQUIRE(r.size() == 3);
  REQUIRE(r == std::vector<double>{0, 0, 1});
  // parallel with different lengths: zero (residual uses unit directions)
  for (double v : c.residuals(lines({2, 0, 0}, {0.1, 0, 0}))) REQUIRE(v == 0.0);
  // antiparallel also satisfies the cross-product residual (documented)
  for (double v : c.residuals(lines({1, 2, 3}, {-2, -4, -6}))) REQUIRE_NEAR(v, 0.0, 1e-15);
  // 30 degrees apart: |r| = sin(30 deg) = 0.5
  r = c.residuals(lines({1, 0, 0}, {std::cos(kPi / 6), std::sin(kPi / 6), 0}));
  REQUIRE_NEAR(std::hypot(r[0], r[1], r[2]), 0.5, 1e-15);
}

TEST_CASE("ParallelConstraint rejects zero-length directions", "[constraint][parallel][errors]") {
  const ParallelConstraint c(kDirA, kDirB);
  // A zero direction would give u x v = 0 and falsely look "parallel".
  REQUIRE_THROWS_AS(c.residuals(lines({0, 0, 0}, {1, 0, 0})), DegenerateGeometryError);
  REQUIRE_THROWS_AS(c.residuals(lines({1, 0, 0}, {0, 0, 0})), DegenerateGeometryError);
  REQUIRE_THROWS_AS(c.residuals(lines({1e-13, 0, 0}, {1, 0, 0})), DegenerateGeometryError);
  REQUIRE_THROWS_AS(c.residuals(lines({std::nan(""), 0, 0}, {1, 0, 0})), DegenerateGeometryError);
}

TEST_CASE("PerpendicularConstraint: r = u . v", "[constraint][perpendicular]") {
  const PerpendicularConstraint c(kDirA, kDirB);
  REQUIRE(c.residualCount() == 1);
  REQUIRE(c.residuals(lines({1, 0, 0}, {0, 7, 0}))[0] == 0.0);
  REQUIRE_NEAR(c.residuals(lines({1, 0, 0}, {1, 1, 0}))[0], std::sqrt(0.5), 1e-15);
  REQUIRE_NEAR(c.residuals(lines({1, 0, 0}, {-3, 0, 0}))[0], -1.0, 1e-15);
  REQUIRE_THROWS_AS(c.residuals(lines({0, 0, 0}, {1, 0, 0})), DegenerateGeometryError);
  REQUIRE_THROWS_AS(c.residuals(lines({1, 0, 0}, {0, 0, 0})), DegenerateGeometryError);
}

TEST_CASE("AngleConstraint: r = u . v - cos(theta)", "[constraint][angle]") {
  const Vector3 u{1, 0, 0};
  const Vector3 at60{5 * std::cos(kPi / 3), 5 * std::sin(kPi / 3), 0};  // 60 degrees from u
  REQUIRE(AngleConstraint(kDirA, kDirB, kPi / 3).residualCount() == 1);
  REQUIRE_NEAR(AngleConstraint(kDirA, kDirB, kPi / 3).residuals(lines(u, at60))[0], 0.0, 1e-15);
  // wrong target: residual is cos(60) - cos(90) = 0.5
  REQUIRE_NEAR(AngleConstraint(kDirA, kDirB, kPi / 2).residuals(lines(u, at60))[0], 0.5, 1e-15);
  REQUIRE_NEAR(AngleConstraint(kDirA, kDirB, 0.0).residuals(lines(u, at60))[0], -0.5, 1e-15);
}

TEST_CASE("AngleConstraint direction convention and antiparallel handling", "[constraint][angle]") {
  const Vector3 u{1, 0, 0};
  const Vector3 same{4, 0, 0};
  const Vector3 opposite{-4, 0, 0};
  // theta = 0: same orientation satisfied, antiparallel is NOT (residual -2).
  REQUIRE_NEAR(AngleConstraint(kDirA, kDirB, 0.0).residuals(lines(u, same))[0], 0.0, 1e-16);
  REQUIRE_NEAR(AngleConstraint(kDirA, kDirB, 0.0).residuals(lines(u, opposite))[0], -2.0, 1e-16);
  // theta = pi: antiparallel satisfied, same orientation is not (residual +2).
  REQUIRE_NEAR(AngleConstraint(kDirA, kDirB, kPi).residuals(lines(u, opposite))[0], 0.0, 1e-16);
  REQUIRE_NEAR(AngleConstraint(kDirA, kDirB, kPi).residuals(lines(u, same))[0], 2.0, 1e-16);
  // Reversing one DirectionRef maps theta -> pi - theta: 60 deg becomes 120 deg.
  const Vector3 at60{std::cos(kPi / 3), std::sin(kPi / 3), 0};
  const AngleConstraint sixty(kDirA, kDirB, kPi / 3);
  const AngleConstraint sixtyReversedB(kDirA, DirectionRef{P3, P2}, kPi / 3);
  REQUIRE_NEAR(sixty.residuals(lines(u, at60))[0], 0.0, 1e-15);
  REQUIRE_NEAR(sixtyReversedB.residuals(lines(u, at60))[0], std::cos(2 * kPi / 3) - 0.5, 1e-15);
  // The angle is unsigned: swapping the two lines does not change it.
  const AngleConstraint swapped(kDirB, kDirA, kPi / 3);
  REQUIRE_NEAR(swapped.residuals(lines(u, at60))[0], 0.0, 1e-15);
}

TEST_CASE("AngleConstraint validation", "[constraint][angle][errors]") {
  REQUIRE_THROWS_AS(AngleConstraint(kDirA, kDirB, -0.1), InvalidArgumentError);
  REQUIRE_THROWS_AS(AngleConstraint(kDirA, kDirB, kPi + 0.1), InvalidArgumentError);
  REQUIRE_THROWS_AS(AngleConstraint(kDirA, kDirB, std::nan("")), InvalidArgumentError);
  REQUIRE_NOTHROW(AngleConstraint(kDirA, kDirB, 0.0));
  REQUIRE_NOTHROW(AngleConstraint(kDirA, kDirB, kPi));
  const AngleConstraint c(kDirA, kDirB, kPi / 4);
  REQUIRE_THROWS_AS(c.residuals(lines({0, 0, 0}, {1, 0, 0})), DegenerateGeometryError);
  REQUIRE_THROWS_AS(c.residuals(lines({1, 0, 0}, {0, 0, 0})), DegenerateGeometryError);
}

TEST_CASE("Constraints validate point ids and buffers", "[constraint][errors]") {
  const std::vector<Point3> pts = {{0, 0, 0}, {1, 0, 0}};
  const DistanceConstraint outOfRange(P0, PointId{7}, 1.0);
  REQUIRE_THROWS_AS(outOfRange.residuals(pts), InvalidArgumentError);

  const CoincidenceConstraint c(P0, P1);
  std::vector<double> tooSmall(2);
  REQUIRE_THROWS_AS(c.evaluate(pts, tooSmall), InvalidArgumentError);
  std::vector<double> tooBig(4);
  REQUIRE_THROWS_AS(c.evaluate(pts, tooBig), InvalidArgumentError);

  REQUIRE(c.referencedPoints() == std::vector<PointId>{P0, P1});
  REQUIRE(ParallelConstraint(kDirA, kDirB).referencedPoints().size() == 4);
}

// ---------------------------------------------------------------------------
// Problem: stable ids, deterministic residual layout, encode/decode
// ---------------------------------------------------------------------------

TEST_CASE("Problem concatenates residuals deterministically", "[problem]") {
  Problem p;
  const PointId a = p.addPoint({0, 0, 0}, "A");
  const PointId b = p.addPoint({3, 4, 0}, "B");
  REQUIRE(a.index == 0);
  REQUIRE(b.index == 1);
  REQUIRE(p.pointName(b) == "B");

  const std::size_t c0 = p.add<FixedPointConstraint>(a, Point3{1, 0, 0});     // 3 residuals
  const std::size_t c1 = p.add<DistanceConstraint>(a, b, 2.0);                // 1 residual
  const std::size_t c2 = p.add<CoincidenceConstraint>(a, b);                  // 3 residuals
  REQUIRE(c0 == 0);
  REQUIRE(c1 == 1);
  REQUIRE(c2 == 2);
  REQUIRE(p.residualOffset(0) == 0);
  REQUIRE(p.residualOffset(1) == 3);
  REQUIRE(p.residualOffset(2) == 4);
  REQUIRE(p.residualCount() == 7);
  REQUIRE(p.variableCount() == 6);

  const Eigen::VectorXd r = p.evaluateResiduals(p.encode(p.initialPoints()));
  REQUIRE(r.size() == 7);
  const std::vector<double> expected = {-1, 0, 0,   // fixed: (0,0,0) - (1,0,0)
                                        3,          // |(3,4,0)| - 2
                                        -3, -4, 0}; // coincident: a - b
  for (Eigen::Index i = 0; i < r.size(); ++i) REQUIRE_NEAR(r[i], expected[static_cast<std::size_t>(i)], 1e-15);
  // repeated evaluation is bit-identical
  REQUIRE(p.evaluateResiduals(p.encode(p.initialPoints())) == r);
}

TEST_CASE("Problem state encoding round-trips", "[problem]") {
  Problem p;
  p.addPoint({1, 2, 3});
  p.addPoint({-4, 5, -6});
  const Eigen::VectorXd x = p.encode(p.initialPoints());
  REQUIRE(x.size() == 6);
  REQUIRE(x[0] == 1);
  REQUIRE(x[3] == -4);
  REQUIRE(x[5] == -6);
  const auto pts = p.decode(x);
  REQUIRE(pts[1] == Point3{-4, 5, -6});
  REQUIRE(p.pointName(PointId{0}) == "p0");  // default name
  REQUIRE_THROWS_AS(p.decode(Eigen::VectorXd(5)), InvalidArgumentError);
}

TEST_CASE("Problem rejects invalid constraints", "[problem][errors]") {
  Problem p;
  const PointId a = p.addPoint({0, 0, 0});
  REQUIRE_THROWS_AS(p.addConstraint(nullptr), InvalidArgumentError);
  // references a point that does not exist (yet)
  REQUIRE_THROWS_AS(p.add<DistanceConstraint>(a, PointId{1}, 1.0), InvalidArgumentError);
  REQUIRE(p.constraintCount() == 0);  // failed adds leave the problem unchanged
  REQUIRE_FALSE(p.validate().empty()); // no constraints
}

// ---------------------------------------------------------------------------
// Finite-difference Jacobian vs. independently derived analytic derivatives
// ---------------------------------------------------------------------------

namespace {

Problem pointsOnly(const std::vector<Point3>& pts) {
  Problem p;
  for (const Point3& q : pts) p.addPoint(q);
  return p;
}

Eigen::Vector3d toEigen(const Vector3& v) { return {v.x, v.y, v.z}; }

Eigen::Matrix3d skew(const Eigen::Vector3d& a) {
  Eigen::Matrix3d s;
  s << 0, -a.z(), a.y(), a.z(), 0, -a.x(), -a.y(), a.x(), 0;
  return s;
}

}  // namespace

TEST_CASE("FD Jacobian of a distance constraint matches (p-q)/|p-q|", "[jacobian]") {
  Problem p = pointsOnly({{1, 2, 3}, {4, -2, 15}});
  p.add<DistanceConstraint>(P0, P1, 1.0);
  const Eigen::VectorXd x = p.encode(p.initialPoints());
  const Eigen::Vector3d g = toEigen((Point3{1, 2, 3} - Point3{4, -2, 15})) / 13.0;  // |d| = 13
  for (JacobianScheme s : {JacobianScheme::Central, JacobianScheme::Forward}) {
    const Eigen::MatrixXd j = finiteDifferenceJacobian(p, x, 1e-6, s);
    REQUIRE(j.rows() == 1);
    REQUIRE(j.cols() == 6);
    const double tol = s == JacobianScheme::Central ? 1e-8 : 1e-5;
    for (int k = 0; k < 3; ++k) {
      REQUIRE_NEAR(j(0, k), g[k], tol);
      REQUIRE_NEAR(j(0, 3 + k), -g[k], tol);
    }
  }
}

TEST_CASE("FD Jacobian of linear constraints is exact", "[jacobian]") {
  Problem p = pointsOnly({{1, 2, 3}, {4, -2, 15}});
  p.add<CoincidenceConstraint>(P0, P1);
  p.add<PointOnPlaneConstraint>(P1, Plane3::fromPointNormal({0, 0, 0}, {0, 3, 4}));
  const Eigen::MatrixXd j = finiteDifferenceJacobian(p, p.encode(p.initialPoints()));
  Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(4, 6);
  expected.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
  expected.block<3, 3>(0, 3) = -Eigen::Matrix3d::Identity();
  expected.block<1, 3>(3, 3) = Eigen::RowVector3d(0, 0.6, 0.8);
  REQUIRE((j - expected).cwiseAbs().maxCoeff() < 1e-8);
}

TEST_CASE("FD Jacobian of an angle constraint matches the analytic gradient", "[jacobian]") {
  const std::vector<Point3> pts = {{0, 0, 0}, {2, 1, 0.5}, {1, 3, -1}, {-1, 5, 2}};
  Problem p = pointsOnly(pts);
  p.add<AngleConstraint>(kDirA, kDirB, 1.0);
  const Eigen::MatrixXd j = finiteDifferenceJacobian(p, p.encode(p.initialPoints()));

  const Eigen::Vector3d u = toEigen(pts[1] - pts[0]);
  const Eigen::Vector3d v = toEigen(pts[3] - pts[2]);
  const Eigen::Vector3d uh = u.normalized();
  const Eigen::Vector3d vh = v.normalized();
  const double c = uh.dot(vh);
  const Eigen::Vector3d dEnd1 = (vh - c * uh) / u.norm();  // d r / d P1 ; P0 gets the negative
  const Eigen::Vector3d dEnd2 = (uh - c * vh) / v.norm();  // d r / d P3 ; P2 gets the negative
  for (int k = 0; k < 3; ++k) {
    REQUIRE_NEAR(j(0, k), -dEnd1[k], 1e-8);
    REQUIRE_NEAR(j(0, 3 + k), dEnd1[k], 1e-8);
    REQUIRE_NEAR(j(0, 6 + k), -dEnd2[k], 1e-8);
    REQUIRE_NEAR(j(0, 9 + k), dEnd2[k], 1e-8);
  }
}

TEST_CASE("FD Jacobian of a parallel constraint matches the analytic gradient", "[jacobian]") {
  const std::vector<Point3> pts = {{0, 0, 0}, {2, 1, 0.5}, {1, 3, -1}, {-1, 5, 2}};
  Problem p = pointsOnly(pts);
  p.add<ParallelConstraint>(kDirA, kDirB);
  const Eigen::MatrixXd j = finiteDifferenceJacobian(p, p.encode(p.initialPoints()));

  const Eigen::Vector3d u = toEigen(pts[1] - pts[0]);
  const Eigen::Vector3d v = toEigen(pts[3] - pts[2]);
  const Eigen::Vector3d uh = u.normalized();
  const Eigen::Vector3d vh = v.normalized();
  const Eigen::Matrix3d pu = (Eigen::Matrix3d::Identity() - uh * uh.transpose()) / u.norm();
  const Eigen::Matrix3d pv = (Eigen::Matrix3d::Identity() - vh * vh.transpose()) / v.norm();
  // r = uh x vh = skew(uh) vh = -skew(vh) uh
  const Eigen::Matrix3d dEnd1 = -skew(vh) * pu;
  const Eigen::Matrix3d dEnd2 = skew(uh) * pv;
  Eigen::MatrixXd expected(3, 12);
  expected << -dEnd1, dEnd1, -dEnd2, dEnd2;
  REQUIRE((j - expected).cwiseAbs().maxCoeff() < 1e-8);
}

TEST_CASE("FD Jacobian validates its step", "[jacobian][errors]") {
  Problem p = pointsOnly({{0, 0, 0}, {1, 0, 0}});
  p.add<DistanceConstraint>(P0, P1, 1.0);
  const Eigen::VectorXd x = p.encode(p.initialPoints());
  REQUIRE_THROWS_AS(finiteDifferenceJacobian(p, x, 0.0), InvalidArgumentError);
  REQUIRE_THROWS_AS(finiteDifferenceJacobian(p, x, -1e-6), InvalidArgumentError);
  REQUIRE_THROWS_AS(finiteDifferenceJacobian(p, x, std::nan("")), InvalidArgumentError);
}
