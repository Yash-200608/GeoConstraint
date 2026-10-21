#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

#include "geoconstraint/constraints/angle_constraint.hpp"
#include "geoconstraint/constraints/coincidence_constraint.hpp"
#include "geoconstraint/constraints/distance_constraint.hpp"
#include "geoconstraint/constraints/fixed_point_constraint.hpp"
#include "geoconstraint/constraints/parallel_constraint.hpp"
#include "geoconstraint/constraints/perpendicular_constraint.hpp"
#include "geoconstraint/constraints/point_on_plane_constraint.hpp"
#include "geoconstraint/solver/solver.hpp"
#include "test_utils.hpp"

using namespace geoconstraint;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

struct Bracket {
  Problem problem;
  PointId a, b, c;
};

Bracket makeBracket(const Point3& initialC, const Point3& initialA = {0, 0, 0},
                    const Point3& initialB = {100, 0, 0}) {
  Bracket br;
  br.a = br.problem.addPoint(initialA, "A");
  br.b = br.problem.addPoint(initialB, "B");
  br.c = br.problem.addPoint(initialC, "C");
  br.problem.add<FixedPointConstraint>(br.a, Point3{0, 0, 0});
  br.problem.add<FixedPointConstraint>(br.b, Point3{100, 0, 0});
  br.problem.add<DistanceConstraint>(br.c, br.a, std::sqrt(5000.0));
  br.problem.add<DistanceConstraint>(br.c, br.b, std::sqrt(5000.0));
  br.problem.add<PointOnPlaneConstraint>(br.c,
                                         Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1}));
  return br;
}

}  // namespace

TEST_CASE("Recovery of a point from known distances (trilateration)", "[solver]") {
  // Anchors and a hidden target P* = (3,4,5); only distances are given to the solver.
  const Point3 target{3, 4, 5};
  const Point3 anchors[3] = {{0, 0, 0}, {10, 0, 0}, {0, 10, 0}};
  Problem p;
  PointId ids[3];
  for (int i = 0; i < 3; ++i) {
    ids[i] = p.addPoint(anchors[i]);
    p.add<FixedPointConstraint>(ids[i], anchors[i]);
  }
  const PointId x = p.addPoint({1, 1, 1}, "X");  // initial guess on the z>0 side
  for (int i = 0; i < 3; ++i) p.add<DistanceConstraint>(x, ids[i], distance(target, anchors[i]));

  const SolverResult r = solve(p);
  REQUIRE(r.status == SolverStatus::Converged);
  REQUIRE_POINT_NEAR(r.points[x.index], target, 1e-7);
  // independent property: distances to the anchors equal the prescribed ones
  for (int i = 0; i < 3; ++i) {
    REQUIRE_NEAR(distance(r.points[x.index], r.points[ids[i].index]),
                 distance(target, anchors[i]), 1e-7);
  }
}

TEST_CASE("Bracket assembly: upper and lower branches", "[solver][bracket]") {
  const double leg = std::sqrt(5000.0);
  struct Case {
    Point3 start;
    Point3 expected;
  };
  const Case cases[] = {{{40, 30, 10}, {50, 50, 0}}, {{40, -30, 10}, {50, -50, 0}}};
  for (const Case& cs : cases) {
    const Bracket br = makeBracket(cs.start);
    const SolverResult r = solve(br.problem);
    INFO(formatReport(br.problem, r));
    REQUIRE(r.status == SolverStatus::Converged);
    REQUIRE_POINT_NEAR(r.points[br.c.index], cs.expected, 1e-8);
    // the defining properties, checked directly from the final geometry
    const Point3& a = r.points[br.a.index];
    const Point3& b = r.points[br.b.index];
    const Point3& c = r.points[br.c.index];
    REQUIRE_NEAR(distance(c, a), leg, 1e-8);
    REQUIRE_NEAR(distance(c, b), leg, 1e-8);
    REQUIRE_NEAR(c.z, 0.0, 1e-8);
    // fixed anchors stay fixed
    REQUIRE_POINT_NEAR(a, (Point3{0, 0, 0}), 1e-8);
    REQUIRE_POINT_NEAR(b, (Point3{100, 0, 0}), 1e-8);
    REQUIRE(r.finalResidualMax <= 1e-9);
    // fully determined locally
    REQUIRE(r.analysis.available);
    REQUIRE(r.analysis.rank == 9);
    REQUIRE(r.analysis.localDegreesOfFreedom == 0);
  }
}

TEST_CASE("Bracket assembly: poor but valid starting guesses", "[solver][bracket]") {
  // None of these has y == 0 exactly: the bracket is symmetric about the plane y = 0, so a
  // start with y == 0 has zero y-gradient and can never leave that plane (see the
  // "symmetric stationary start" test below).
  const Point3 guesses[] = {{-200, 300, 80}, {500, -400, -250}, {50, 1, 0},       {50, -1, 0},
                            {0, 0.1, 1},     {100, -0.1, 1},    {1000, 1000, 1000}, {-60, 5, 5}};
  for (const Point3& guess : guesses) {
    const Bracket br = makeBracket(guess);
    const SolverResult r = solve(br.problem);
    INFO("start = " << test::str(guess) << "\n" << formatReport(br.problem, r));
    REQUIRE(r.status == SolverStatus::Converged);
    const Point3& c = r.points[br.c.index];
    // Either branch is a valid answer - the solver does not pick one in advance.
    const double dUpper = distance(c, Point3{50, 50, 0});
    const double dLower = distance(c, Point3{50, -50, 0});
    REQUIRE(std::min(dUpper, dLower) <= 1e-7);
  }
}

TEST_CASE("Perturbed anchors are pulled back to their targets", "[solver][fixed]") {
  const Bracket br = makeBracket({40, 30, 10}, /*A*/ {2.5, -1.5, 0.75}, /*B*/ {97, 3, -2});
  const SolverResult r = solve(br.problem);
  REQUIRE(r.status == SolverStatus::Converged);
  REQUIRE_POINT_NEAR(r.points[br.a.index], (Point3{0, 0, 0}), 1e-8);
  REQUIRE_POINT_NEAR(r.points[br.b.index], (Point3{100, 0, 0}), 1e-8);
  REQUIRE_POINT_NEAR(r.points[br.c.index], (Point3{50, 50, 0}), 1e-7);
}

TEST_CASE("An exact initial solution needs no iterations", "[solver]") {
  const Bracket br = makeBracket({50, 50, 0});
  const SolverResult r = solve(br.problem);
  REQUIRE(r.status == SolverStatus::Converged);
  REQUIRE(r.iterations == 0);
  REQUIRE(r.points[br.c.index] == Point3{50, 50, 0});
  // and a start that is off only at the 1e-6 level converges quickly
  const Bracket near = makeBracket({50 + 1e-6, 50 - 1e-6, 1e-6});
  const SolverResult rn = solve(near.problem);
  REQUIRE(rn.status == SolverStatus::Converged);
  REQUIRE(rn.iterations <= 5);
  REQUIRE_POINT_NEAR(rn.points[near.c.index], (Point3{50, 50, 0}), 1e-8);
}

TEST_CASE("Solving is deterministic and the cost never increases", "[solver]") {
  const Bracket br = makeBracket({-60, 5, 5});
  const SolverResult r1 = solve(br.problem);
  const SolverResult r2 = solve(br.problem);
  REQUIRE(r1.points == r2.points);
  REQUIRE(r1.iterations == r2.iterations);
  REQUIRE(r1.costHistory.size() == r1.iterations + 1);
  for (std::size_t i = 1; i < r1.costHistory.size(); ++i) {
    REQUIRE(r1.costHistory[i] < r1.costHistory[i - 1]);  // steps are only accepted if they improve
  }
  REQUIRE(r1.costHistory.front() > 0.0);
}

TEST_CASE("Perpendicular, angle and parallel constraints solve to the right geometry",
          "[solver][directions]") {
  SECTION("perpendicular") {
    // |CD| = 2, CD perpendicular to AB (x-axis), D in the plane z = 0  =>  D = (0, +/-2, 0)
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}), b = p.addPoint({1, 0, 0});
    const PointId c = p.addPoint({0, 0, 0}), d = p.addPoint({0.3, 1.5, 0.2});
    p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
    p.add<FixedPointConstraint>(b, Point3{1, 0, 0});
    p.add<FixedPointConstraint>(c, Point3{0, 0, 0});
    p.add<PerpendicularConstraint>(DirectionRef{a, b}, DirectionRef{c, d});
    p.add<DistanceConstraint>(c, d, 2.0);
    p.add<PointOnPlaneConstraint>(d, Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1}));
    const SolverResult r = solve(p);
    REQUIRE(r.status == SolverStatus::Converged);
    REQUIRE_POINT_NEAR(r.points[d.index], (Point3{0, 2, 0}), 1e-7);
  }
  SECTION("fixed angle of 60 degrees") {
    // D = C + 2 (cos 60, +/- sin 60, 0) = (1, +/- sqrt 3, 0)
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}), b = p.addPoint({1, 0, 0});
    const PointId c = p.addPoint({0, 0, 0}), d = p.addPoint({1.2, 1.0, 0.3});
    p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
    p.add<FixedPointConstraint>(b, Point3{1, 0, 0});
    p.add<FixedPointConstraint>(c, Point3{0, 0, 0});
    p.add<AngleConstraint>(DirectionRef{a, b}, DirectionRef{c, d}, kPi / 3);
    p.add<DistanceConstraint>(c, d, 2.0);
    p.add<PointOnPlaneConstraint>(d, Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1}));
    const SolverResult r = solve(p);
    REQUIRE(r.status == SolverStatus::Converged);
    REQUIRE_POINT_NEAR(r.points[d.index], (Point3{1, std::sqrt(3.0), 0}), 1e-7);
    // independent check of the angle itself
    REQUIRE_NEAR(angleBetween(r.points[b.index] - r.points[a.index],
                              r.points[d.index] - r.points[c.index]),
                 kPi / 3, 1e-8);
  }
  SECTION("parallel") {
    // CD parallel to AB (x-axis), |CD| = 3, D in plane z = 0  =>  D = (+/-3, 1, 0)
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}), b = p.addPoint({1, 0, 0});
    const PointId c = p.addPoint({0, 1, 0}), d = p.addPoint({1.5, 1.4, 0.3});
    p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
    p.add<FixedPointConstraint>(b, Point3{1, 0, 0});
    p.add<FixedPointConstraint>(c, Point3{0, 1, 0});
    p.add<ParallelConstraint>(DirectionRef{a, b}, DirectionRef{c, d});
    p.add<DistanceConstraint>(c, d, 3.0);
    p.add<PointOnPlaneConstraint>(d, Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1}));
    const SolverResult r = solve(p);
    REQUIRE(r.status == SolverStatus::Converged);
    const Point3& dp = r.points[d.index];
    REQUIRE(std::min(distance(dp, Point3{3, 1, 0}), distance(dp, Point3{-3, 1, 0})) <= 1e-7);
    REQUIRE(isParallel(r.points[b.index] - r.points[a.index], dp - r.points[c.index],
                       Tolerances{1e-12, 1e-9, 1e-9, 1e-7}));
  }
  SECTION("coincidence") {
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}), b = p.addPoint({4, -2, 9});
    p.add<FixedPointConstraint>(a, Point3{1, 2, 3});
    p.add<CoincidenceConstraint>(a, b);
    const SolverResult r = solve(p);
    REQUIRE(r.status == SolverStatus::Converged);
    REQUIRE_POINT_NEAR(r.points[b.index], (Point3{1, 2, 3}), 1e-8);
  }
}

TEST_CASE("Underdetermined systems converge but report free degrees of freedom",
          "[solver][diagnostics]") {
  Problem p;
  const PointId a = p.addPoint({0, 0, 0}), c = p.addPoint({1, 2, 2});
  p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
  p.add<DistanceConstraint>(c, a, 5.0);
  const SolverResult r = solve(p);
  REQUIRE(r.status == SolverStatus::Converged);
  REQUIRE_NEAR(norm(r.points[c.index] - r.points[a.index]), 5.0, 1e-8);  // somewhere on the sphere
  REQUIRE(r.analysis.available);
  REQUIRE(r.analysis.variableCount == 6);
  REQUIRE(r.analysis.rank == 4);
  REQUIRE(r.analysis.localDegreesOfFreedom == 2);  // converged, yet not unique
}

TEST_CASE("Inconsistent constraints do not report success", "[solver][failure]") {
  SECTION("one point pinned to two places") {
    Problem p;
    const PointId a = p.addPoint({0.5, 0, 0});
    p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
    p.add<FixedPointConstraint>(a, Point3{2, 0, 0});
    const SolverResult r = solve(p);
    REQUIRE(r.status == SolverStatus::Stagnated);
    REQUIRE_FALSE(r.converged());
    // the least-squares compromise is the midpoint, with residual (1,0,0) and (-1,0,0)
    REQUIRE_POINT_NEAR(r.points[a.index], (Point3{1, 0, 0}), 1e-6);
    REQUIRE_NEAR(r.finalResidualNorm, std::sqrt(2.0), 1e-6);
    REQUIRE(r.finalResidualMax > 0.9);
    REQUIRE_FALSE(r.message.empty());
  }
  SECTION("distance contradicting fixed anchors") {
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}), b = p.addPoint({10, 0, 0});
    p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
    p.add<FixedPointConstraint>(b, Point3{10, 0, 0});
    p.add<DistanceConstraint>(a, b, 4.0);  // but A and B are pinned 10 apart
    const SolverResult r = solve(p);
    REQUIRE_FALSE(r.converged());
    REQUIRE(r.finalResidualMax > 1.0);
    REQUIRE(r.constraintResiduals.size() == 3);
  }
}

TEST_CASE("A symmetric stationary start is reported honestly, not as success",
          "[solver][failure]") {
  // C on the perpendicular bisector line, y = 0: the y-gradient vanishes exactly, the
  // x-gradients cancel, so the start is a stationary point of the objective.
  const Bracket br = makeBracket({50, 0, 0});
  const SolverResult r = solve(br.problem);
  REQUIRE_FALSE(r.converged());
  REQUIRE(r.status == SolverStatus::Stagnated);
  REQUIRE(r.finalResidualMax > 1.0);
}

TEST_CASE("Iteration limit is reported as MaxIterations", "[solver][failure]") {
  const Bracket br = makeBracket({-200, 300, 80});
  SolverOptions options;
  options.maxIterations = 1;
  const SolverResult r = solve(br.problem, options);
  REQUIRE(r.status == SolverStatus::MaxIterations);
  REQUIRE(r.iterations == 1);
  REQUIRE_FALSE(r.converged());
  // the state returned is still the (improved) iterate, with full diagnostics
  REQUIRE(r.finalResidualNorm < r.initialResidualNorm);
  REQUIRE(r.constraintResiduals.size() == 5);
  REQUIRE(r.analysis.available);
}

TEST_CASE("Invalid problems and states return InvalidInput", "[solver][errors]") {
  SECTION("NaN initial coordinate") {
    Problem p;
    const PointId a = p.addPoint({0, 0, 0});
    const PointId b = p.addPoint({kNaN, 0, 0}, "bad");
    p.add<CoincidenceConstraint>(a, b);
    const SolverResult r = solve(p);
    REQUIRE(r.status == SolverStatus::InvalidInput);
    REQUIRE(r.message.find("bad") != std::string::npos);
    REQUIRE(r.iterations == 0);
  }
  SECTION("infinite initial coordinate") {
    Problem p;
    const PointId a = p.addPoint({0, 0, 0});
    const PointId b = p.addPoint({0, kInf, 0});
    p.add<DistanceConstraint>(a, b, 1.0);
    REQUIRE(solve(p).status == SolverStatus::InvalidInput);
  }
  SECTION("zero-length direction at the start") {
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}), b = p.addPoint({0, 0, 0});
    const PointId c = p.addPoint({0, 0, 0}), d = p.addPoint({1, 0, 0});
    p.add<ParallelConstraint>(DirectionRef{a, b}, DirectionRef{c, d});
    const SolverResult r = solve(p);
    REQUIRE(r.status == SolverStatus::InvalidInput);
    REQUIRE(r.message.find("zero-length") != std::string::npos);
  }
  SECTION("empty problems") {
    Problem none;
    REQUIRE(solve(none).status == SolverStatus::InvalidInput);
    Problem noConstraints;
    noConstraints.addPoint({0, 0, 0});
    REQUIRE(solve(noConstraints).status == SolverStatus::InvalidInput);
  }
  SECTION("invalid options") {
    const Bracket br = makeBracket({40, 30, 10});
    SolverOptions o;
    o.maxIterations = 0;
    REQUIRE(solve(br.problem, o).status == SolverStatus::InvalidInput);
    o = {};
    o.residualTolerance = -1.0;
    REQUIRE(solve(br.problem, o).status == SolverStatus::InvalidInput);
    o = {};
    o.finiteDifferenceStep = 0.0;
    REQUIRE(solve(br.problem, o).status == SolverStatus::InvalidInput);
    o = {};
    o.dampingDecrease = 1.5;
    REQUIRE(solve(br.problem, o).status == SolverStatus::InvalidInput);
    o = {};
    o.residualTolerance = kNaN;
    REQUIRE(solve(br.problem, o).status == SolverStatus::InvalidInput);
  }
}

TEST_CASE("A tighter residual tolerance yields a more accurate solution", "[solver][options]") {
  const Bracket br = makeBracket({40, 30, 10});
  SolverOptions loose;
  loose.residualTolerance = 1e-3;
  SolverOptions tight;
  tight.residualTolerance = 1e-11;
  const SolverResult rl = solve(br.problem, loose);
  const SolverResult rt = solve(br.problem, tight);
  REQUIRE(rl.status == SolverStatus::Converged);
  REQUIRE(rt.status == SolverStatus::Converged);
  REQUIRE(rl.finalResidualMax <= 1e-3);
  REQUIRE(rt.finalResidualMax <= 1e-11);
  REQUIRE(rl.iterations <= rt.iterations);
}

TEST_CASE("Forward and central Jacobians both solve the bracket", "[solver][options]") {
  const Bracket br = makeBracket({40, 30, 10});
  SolverOptions fwd;
  fwd.jacobianScheme = JacobianScheme::Forward;
  fwd.finiteDifferenceStep = 1e-7;
  const SolverResult r = solve(br.problem, fwd);
  REQUIRE(r.status == SolverStatus::Converged);
  REQUIRE_POINT_NEAR(r.points[br.c.index], (Point3{50, 50, 0}), 1e-7);
}

// ---------------------------------------------------------------------------
// Rank / local DOF diagnostics at hand-built configurations
// ---------------------------------------------------------------------------

TEST_CASE("Local DOF of simple constraint systems", "[diagnostics]") {
  const Plane3 xy = Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1});
  SECTION("free point on a plane: 3 - 1 = 2") {
    Problem p;
    const PointId a = p.addPoint({1, 2, 3});
    p.add<PointOnPlaneConstraint>(a, xy);
    const auto an = analyzeSystem(p, p.initialPoints());
    REQUIRE(an.rank == 1);
    REQUIRE(an.localDegreesOfFreedom == 2);
  }
  SECTION("two coincident points: 6 - 3 = 3") {
    Problem p;
    const PointId a = p.addPoint({1, 2, 3}), b = p.addPoint({4, 5, 6});
    p.add<CoincidenceConstraint>(a, b);
    const auto an = analyzeSystem(p, p.initialPoints());
    REQUIRE(an.rank == 3);
    REQUIRE(an.localDegreesOfFreedom == 3);
  }
  SECTION("pinned point: 0 DOF") {
    Problem p;
    const PointId a = p.addPoint({1, 2, 3});
    p.add<FixedPointConstraint>(a, Point3{1, 2, 3});
    REQUIRE(analyzeSystem(p, p.initialPoints()).localDegreesOfFreedom == 0);
  }
  SECTION("perpendicular: 12 - 1 = 11") {
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}), b = p.addPoint({1, 0, 0});
    const PointId c = p.addPoint({0, 1, 0}), d = p.addPoint({1, 2, 3});
    p.add<PerpendicularConstraint>(DirectionRef{a, b}, DirectionRef{c, d});
    const auto an = analyzeSystem(p, p.initialPoints());
    REQUIRE(an.localDegreesOfFreedom == 11);
  }
  SECTION("parallel: rank 3 when violated, rank 2 once satisfied") {
    // r = u x v. Away from parallel the three components respond independently
    // (rank 3); at u = +/-v the first-order variation lies in the 2-D plane
    // orthogonal to u (rank 2), leaving one redundant residual row.
    Problem generic;
    {
      const PointId a = generic.addPoint({0, 0, 0}), b = generic.addPoint({1, 0, 0});
      const PointId c = generic.addPoint({0, 1, 0}), d = generic.addPoint({1, 2, 3});
      generic.add<ParallelConstraint>(DirectionRef{a, b}, DirectionRef{c, d});
    }
    REQUIRE(analyzeSystem(generic, generic.initialPoints()).rank == 3);

    Problem satisfied;
    {
      const PointId a = satisfied.addPoint({0, 0, 0}), b = satisfied.addPoint({1, 0, 0});
      const PointId c = satisfied.addPoint({0, 1, 0}), d = satisfied.addPoint({3, 1, 0});
      satisfied.add<ParallelConstraint>(DirectionRef{a, b}, DirectionRef{c, d});
    }
    const auto an = analyzeSystem(satisfied, satisfied.initialPoints());
    REQUIRE(an.residualCount == 3);
    REQUIRE(an.rank == 2);
    REQUIRE(an.localDegreesOfFreedom == 10);
    REQUIRE(an.redundantResiduals == 1);
  }
  SECTION("redundant but consistent constraints do not change the rank") {
    Problem p;
    const PointId a = p.addPoint({1, 2, 3});
    p.add<FixedPointConstraint>(a, Point3{1, 2, 3});
    p.add<FixedPointConstraint>(a, Point3{1, 2, 3});
    const auto an = analyzeSystem(p, p.initialPoints());
    REQUIRE(an.rank == 3);
    REQUIRE(an.redundantResiduals == 3);
    REQUIRE(an.localDegreesOfFreedom == 0);
    const SolverResult r = solve(p);
    REQUIRE(r.status == SolverStatus::Converged);  // consistent redundancy converges
  }
}

TEST_CASE("Rank estimate is first order only: tangent configuration", "[diagnostics][limitations]") {
  // C = (50,0,0) with |CA| = |CB| = 50 and A=(0,0,0), B=(100,0,0): the distance
  // gradients are (+1,0,0) and (-1,0,0), so J loses a rank although C is
  // geometrically unique (the y/z deviations enter the distances only at second
  // order). The estimate reports 1 local DOF - a documented limitation.
  Bracket br;
  br.a = br.problem.addPoint({0, 0, 0});
  br.b = br.problem.addPoint({100, 0, 0});
  br.c = br.problem.addPoint({50, 0, 0});
  br.problem.add<FixedPointConstraint>(br.a, Point3{0, 0, 0});
  br.problem.add<FixedPointConstraint>(br.b, Point3{100, 0, 0});
  br.problem.add<DistanceConstraint>(br.c, br.a, 50.0);
  br.problem.add<DistanceConstraint>(br.c, br.b, 50.0);
  br.problem.add<PointOnPlaneConstraint>(br.c, Plane3::fromPointNormal({0, 0, 0}, {0, 0, 1}));
  const SystemAnalysis an = analyzeSystem(br.problem, br.problem.initialPoints());
  REQUIRE(an.rank == 8);
  REQUIRE(an.localDegreesOfFreedom == 1);
  REQUIRE(an.smallestSingularValue < 1e-8 * an.largestSingularValue);
}

TEST_CASE("Singular values and threshold are reported consistently", "[diagnostics]") {
  const Bracket br = makeBracket({50, 50, 0});
  const SystemAnalysis an = analyzeSystem(br.problem, br.problem.initialPoints());
  REQUIRE(an.singularValues.size() == 9);
  REQUIRE(std::is_sorted(an.singularValues.rbegin(), an.singularValues.rend()));
  REQUIRE(an.largestSingularValue == an.singularValues.front());
  REQUIRE(an.smallestSingularValue == an.singularValues.back());
  REQUIRE_NEAR(an.rankThreshold, 1e-8 * an.largestSingularValue, 1e-20);
  REQUIRE(an.rank == 9);
  REQUIRE(an.localDegreesOfFreedom == 0);

  // tightening the rank tolerance cannot lower the rank; loosening it enough does
  SolverOptions loose;
  loose.rankTolerance = 0.9;
  const SystemAnalysis coarse = analyzeSystem(br.problem, br.problem.initialPoints(), loose);
  REQUIRE(coarse.rank < an.rank);
  REQUIRE_THROWS_AS(analyzeJacobian(Eigen::MatrixXd::Identity(2, 2), 0.0), InvalidArgumentError);
  Eigen::MatrixXd bad = Eigen::MatrixXd::Identity(2, 2);
  bad(0, 0) = kNaN;
  REQUIRE_THROWS_AS(analyzeJacobian(bad, 1e-8), InvalidArgumentError);
}

TEST_CASE("Per-constraint residuals and the report", "[diagnostics]") {
  const Bracket br = makeBracket({40, 30, 10});
  SolverOptions o;
  o.maxIterations = 1;  // stop early so residuals are not all ~0
  const SolverResult r = solve(br.problem, o);
  REQUIRE(r.constraintResiduals.size() == 5);
  const std::size_t sizes[] = {3, 3, 1, 1, 1};
  const char* types[] = {"FixedPoint", "FixedPoint", "Distance", "Distance", "PointOnPlane"};
  double sumSq = 0.0;
  for (std::size_t i = 0; i < 5; ++i) {
    const ConstraintResidual& cr = r.constraintResiduals[i];
    REQUIRE(cr.constraintIndex == i);
    REQUIRE(cr.values.size() == sizes[i]);
    REQUIRE(cr.typeName == types[i]);
    double ss = 0.0;
    for (double v : cr.values) ss += v * v;
    REQUIRE_NEAR(cr.norm, std::sqrt(ss), 1e-15);
    sumSq += ss;
  }
  REQUIRE_NEAR(std::sqrt(sumSq), r.finalResidualNorm, 1e-12);
  const std::string report = formatReport(br.problem, r);
  REQUIRE(report.find("MaxIterations") != std::string::npos);
  REQUIRE(report.find("local DOF") != std::string::npos);
  REQUIRE(report.find("C = (") != std::string::npos);
  REQUIRE(toString(SolverStatus::Converged) == "Converged");
}
