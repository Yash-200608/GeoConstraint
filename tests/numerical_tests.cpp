#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

#include "geoconstraint/constraints/distance_constraint.hpp"
#include "geoconstraint/constraints/fixed_point_constraint.hpp"
#include "geoconstraint/constraints/parallel_constraint.hpp"
#include "geoconstraint/constraints/point_on_plane_constraint.hpp"
#include "geoconstraint/geometry/distances.hpp"
#include "geoconstraint/geometry/intersections.hpp"
#include "geoconstraint/geometry/projections.hpp"
#include "geoconstraint/solver/solver.hpp"
#include "test_utils.hpp"

using namespace geoconstraint;

TEST_CASE("approxEqual boundaries follow the documented policy", "[numerical][tolerances]") {
  // |a-b| <= length + relative * max(|a|,|b|) with the defaults 1e-9 and 1e-9.
  REQUIRE(approxEqual(0.0, 0.9e-9));
  REQUIRE_FALSE(approxEqual(0.0, 1.1e-9));
  REQUIRE(approxEqual(1e6, 1e6 + 1.0e-3));       // tolerance ~ 1e-9 + 1e-3
  REQUIRE_FALSE(approxEqual(1e6, 1e6 + 2.0e-3));
  REQUIRE(approxEqual(1e-6, 1e-6 + 0.5e-9));
  REQUIRE_FALSE(approxEqual(1e-6, 1e-6 + 5e-9));  // relative scale does not hide it at tiny sizes
  const double inf = std::numeric_limits<double>::infinity();
  REQUIRE_FALSE(approxEqual(inf, inf));           // non-finite never "equal"
  REQUIRE_FALSE(approxEqual(std::nan(""), 0.0));
  REQUIRE(lengthTolerance(kDefaultTolerances, 0.0) == 1e-9);
  REQUIRE_NEAR(lengthTolerance(kDefaultTolerances, -1e3), 1e-9 + 1e-6, 1e-18);  // sign-independent
}

TEST_CASE("Almost-zero directions are rejected at the documented threshold", "[numerical]") {
  REQUIRE_THROWS_AS(normalized({1e-13, 0, 0}), DegenerateGeometryError);
  REQUIRE_THROWS_AS(normalized({5e-13, 5e-13, 5e-13}), DegenerateGeometryError);  // |v| ~ 8.7e-13
  const Vector3 small = normalized({1e-11, 0, 0});  // above zeroLength: a valid direction
  REQUIRE_NEAR(small.x, 1.0, 1e-15);
  REQUIRE_THROWS_AS(Line3::fromTwoPoints({1, 1, 1}, {1, 1, 1 + 1e-13}), DegenerateGeometryError);
  // a stricter policy can reclassify the same vector
  Tolerances strict;
  strict.zeroLength = 1e-9;
  REQUIRE_THROWS_AS(normalized({1e-11, 0, 0}, strict), DegenerateGeometryError);
}

TEST_CASE("Parallel/perpendicular classification at the angular tolerance", "[numerical]") {
  // For a = x-axis and b = (1, s, 0): sin(angle) ~ s.
  REQUIRE(isParallel({1, 0, 0}, {1, 0.5e-9, 0}));
  REQUIRE_FALSE(isParallel({1, 0, 0}, {1, 2e-9, 0}));
  REQUIRE_FALSE(isParallel({1, 0, 0}, {1, 1e-6, 0}));
  // cos(angle) ~ s for b = (s, 1, 0)
  REQUIRE(isPerpendicular({1, 0, 0}, {0.5e-9, 1, 0}));
  REQUIRE_FALSE(isPerpendicular({1, 0, 0}, {2e-9, 1, 0}));
  // the test is scale-invariant (uses |a||b|)
  REQUIRE(isParallel({1e-6, 0, 0}, {1e6, 0.5e-3, 0}));        // sin = 5e-10
  REQUIRE_FALSE(isParallel({1e-6, 0, 0}, {1e6, 2e-3, 0}));     // sin = 2e-9
}

TEST_CASE("Small angles and near-pi angles stay accurate", "[numerical]") {
  const double eps = 1e-8;
  const double small = angleBetween({1, 0, 0}, {std::cos(eps), std::sin(eps), 0});
  REQUIRE_THAT(small, Catch::Matchers::WithinRel(eps, 1e-7));  // acos would lose ~8 digits here
  const double nearPi = angleBetween({1, 0, 0}, {-std::cos(eps), std::sin(eps), 0});
  REQUIRE_THAT(std::numbers::pi - nearPi, Catch::Matchers::WithinRel(eps, 1e-7));
}

TEST_CASE("Geometry at very small coordinates", "[numerical][scale]") {
  const double s = 1e-6;
  const Line3 line = Line3::fromPointDirection({0, 0, 0}, {1, 0, 0});
  REQUIRE_THAT(distancePointToLine({3 * s, 4 * s, 0}, line),
               Catch::Matchers::WithinRel(4 * s, 1e-12));
  const Plane3 plane = Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1});
  const auto hit = intersectLinePlane(Line3::fromTwoPoints({0, 0, 2 * s}, {s, s, s}), plane);
  REQUIRE(hit.relation == LinePlaneRelation::UniquePoint);
  REQUIRE_THAT(hit.point->x, Catch::Matchers::WithinRel(2 * s, 1e-12));
  // A directed segment of length 1e-6 still defines a direction (>> zeroLength).
  REQUIRE_NOTHROW(Line3::fromTwoPoints({0, 0, 0}, {s, 0, 0}));
}

TEST_CASE("Geometry at large coordinates", "[numerical][scale]") {
  const double big = 1e6;
  const Plane3 plane = Plane3::fromPointNormal({big, big, big}, {0, 0, 1});
  REQUIRE_NEAR(distancePointToPlane({big + 0.123, big, big + 0.5}, plane), 0.5, 1e-9);
  const Line3 line = Line3::fromPointDirection({big, 0, big}, {3, 4, 0});
  // w = (4,3,0), d = (0.6,0.8,0): |w x d| = |4*0.8 - 3*0.6| = 1.4
  REQUIRE_NEAR(distancePointToLine({big + 4, 3, big}, line), 1.4, 1e-9);
  // Projection stays on the line to within the resolution of the coordinates.
  const Point3 proj = projectPointOnLine({big + 10, 7, big + 2}, line);
  REQUIRE_NEAR(distancePointToLine(proj, line), 0.0, 1e-8);
  // Coincident-line classification uses the scale-aware tolerance: 1e-3 off at
  // 1e7 is "the same line", 1e-3 off at the origin is not.
  const Line3 farA = Line3::fromPointDirection({1e7, 0, 0}, {1, 0, 0});
  const Line3 farB = Line3::fromPointDirection({1e7, 1e-3, 0}, {1, 0, 0});
  REQUIRE(closestPointsBetweenLines(farA, farB).relation == LineLineRelation::Coincident);
  const Line3 nearA = Line3::fromPointDirection({0, 0, 0}, {1, 0, 0});
  const Line3 nearB = Line3::fromPointDirection({0, 1e-3, 0}, {1, 0, 0});
  REQUIRE(closestPointsBetweenLines(nearA, nearB).relation == LineLineRelation::Parallel);
}

TEST_CASE("Almost-zero directions in constraints do not satisfy them", "[numerical][constraint]") {
  Problem p;
  const PointId a = p.addPoint({0, 0, 0}), b = p.addPoint({1e-13, 0, 0});  // direction below zeroLength
  const PointId c = p.addPoint({0, 1, 0}), d = p.addPoint({1, 1, 0});
  p.add<ParallelConstraint>(DirectionRef{a, b}, DirectionRef{c, d});
  REQUIRE(solve(p).status == SolverStatus::InvalidInput);

  // A direction just above the threshold is accepted and evaluated normally.
  Problem q;
  const PointId a2 = q.addPoint({0, 0, 0}), b2 = q.addPoint({1e-9, 0, 0});
  const PointId c2 = q.addPoint({0, 1, 0}), d2 = q.addPoint({1, 1, 1});
  q.add<ParallelConstraint>(DirectionRef{a2, b2}, DirectionRef{c2, d2});
  const Eigen::VectorXd r = q.evaluateResiduals(q.encode(q.initialPoints()));
  REQUIRE_NEAR(r.norm(), std::sqrt(0.5), 1e-12);  // x-axis vs (1,0,1)/sqrt2: sin = 1/sqrt2
}

TEST_CASE("Near-parallel directions give residual equal to sin(angle)", "[numerical][constraint]") {
  const double eps = 1e-7;
  Problem p;
  const PointId a = p.addPoint({0, 0, 0}), b = p.addPoint({1, 0, 0});
  const PointId c = p.addPoint({0, 1, 0}), d = p.addPoint({std::cos(eps), 1 + std::sin(eps), 0});
  p.add<ParallelConstraint>(DirectionRef{a, b}, DirectionRef{c, d});
  const Eigen::VectorXd r = p.evaluateResiduals(p.encode(p.initialPoints()));
  REQUIRE_THAT(r.norm(), Catch::Matchers::WithinRel(std::sin(eps), 1e-9));
}

namespace {

// Bracket scaled by `scale` (all lengths multiplied).
SolverResult solveScaledBracket(double scale, double residualTolerance, Point3* cOut) {
  Problem p;
  const PointId a = p.addPoint({0, 0, 0});
  const PointId b = p.addPoint({100 * scale, 0, 0});
  const PointId c = p.addPoint({40 * scale, 30 * scale, 10 * scale});
  p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
  p.add<FixedPointConstraint>(b, Point3{100 * scale, 0, 0});
  p.add<DistanceConstraint>(c, a, std::sqrt(5000.0) * scale);
  p.add<DistanceConstraint>(c, b, std::sqrt(5000.0) * scale);
  p.add<PointOnPlaneConstraint>(c, Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1}));
  SolverOptions options;
  options.residualTolerance = residualTolerance;
  SolverResult r = solve(p, options);
  *cOut = r.points[c.index];
  return r;
}

}  // namespace

TEST_CASE("The solver works across coordinate scales", "[numerical][scale][solver]") {
  for (const double scale : {1e-3, 1.0, 1e3, 1e6}) {
    Point3 c;
    const SolverResult r = solveScaledBracket(scale, 1e-9, &c);
    INFO("scale = " << scale);
    REQUIRE(r.status == SolverStatus::Converged);
    REQUIRE_POINT_NEAR(c, (Point3{50 * scale, 50 * scale, 0}), 1e-7 + 1e-9 * scale);
  }
}

TEST_CASE("Converged always implies the residual criterion holds", "[numerical][solver]") {
  // The tolerance is absolute. Tolerances below what double precision can deliver for
  // the problem must end as a non-success status, never as a false Converged.
  for (const double tol : {1e-6, 1e-9, 1e-12, 1e-15, 1e-30, 1e-300}) {
    SolverOptions options;
    options.residualTolerance = tol;
    Problem p;
    const PointId a = p.addPoint({0, 0, 0});
    const PointId b = p.addPoint({10, 0, 0});
    const PointId x = p.addPoint({1, 1, 1});
    p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
    p.add<FixedPointConstraint>(b, Point3{10, 0, 0});
    p.add<DistanceConstraint>(x, a, std::sqrt(2.0) * 7.0);  // irrational targets
    p.add<DistanceConstraint>(x, b, std::sqrt(3.0) * 5.0);
    p.add<PointOnPlaneConstraint>(x, Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1}));
    const SolverResult r = solve(p, options);
    INFO("tolerance = " << tol << ", status = " << toString(r.status) << ": " << r.message);
    if (r.converged()) REQUIRE(r.finalResidualMax <= tol);
    REQUIRE(r.status != SolverStatus::InvalidInput);
    REQUIRE(r.status != SolverStatus::NumericalFailure);
  }
}

TEST_CASE("Residual tolerance boundary", "[numerical][solver]") {
  auto iterationsFor = [](double zOffset) {
    Problem p;
    const PointId a = p.addPoint({0, 0, 0});
    const PointId c = p.addPoint({3, 4, zOffset});
    p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
    p.add<PointOnPlaneConstraint>(c, Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1}));
    return solve(p);
  };
  const SolverResult inside = iterationsFor(0.5e-9);   // below residualTolerance = 1e-9
  REQUIRE(inside.status == SolverStatus::Converged);
  REQUIRE(inside.iterations == 0);
  const SolverResult outside = iterationsFor(2e-9);    // above the tolerance: must iterate
  REQUIRE(outside.status == SolverStatus::Converged);
  REQUIRE(outside.iterations >= 1);
  REQUIRE(outside.finalResidualMax <= 1e-9);
}

TEST_CASE("Plane through nearly collinear points is rejected by the angular tolerance",
          "[numerical]") {
  // sin(angle between edges) = 5e-13 < 1e-9 -> degenerate
  REQUIRE_THROWS_AS(Plane3::fromThreePoints({0, 0, 0}, {1, 0, 0}, {2, 1e-12, 0}),
                    DegenerateGeometryError);
  // sin = 5e-7 -> accepted, with a well-defined (if sensitive) normal
  const Plane3 p = Plane3::fromThreePoints({0, 0, 0}, {1, 0, 0}, {2, 1e-6, 0});
  REQUIRE_NEAR(std::abs(p.normal().z), 1.0, 1e-12);
}
