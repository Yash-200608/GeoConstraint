// Bracket assembly: the primary end-to-end demonstration.
//
//   A = (0,0,0) and B = (100,0,0) are anchored with fixed-point constraints.
//   C is free, constrained by |C-A| = |C-B| = sqrt(5000) and z_C = 0.
//
// The two valid solutions are C = (50, +50, 0) and C = (50, -50, 0). Which one
// the solver finds depends only on the starting guess; nothing is hard-coded.
//
// Usage: bracket_assembly            runs both documented starting guesses
//        bracket_assembly x y z      runs a single custom starting guess for C

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include "geoconstraint/constraints/distance_constraint.hpp"
#include "geoconstraint/constraints/fixed_point_constraint.hpp"
#include "geoconstraint/constraints/point_on_plane_constraint.hpp"
#include "geoconstraint/solver/solver.hpp"

using namespace geoconstraint;

namespace {

constexpr double kCheckTolerance = 1e-8;  // acceptance threshold for the self-check

struct Bracket {
  Problem problem;
  PointId a, b, c;
};

Bracket buildBracket(const Point3& initialC) {
  Bracket br;
  const Point3 a0{0, 0, 0};
  const Point3 b0{100, 0, 0};
  br.a = br.problem.addPoint(a0, "A");
  br.b = br.problem.addPoint(b0, "B");
  br.c = br.problem.addPoint(initialC, "C");

  const double legLength = std::sqrt(5000.0);
  br.problem.add<FixedPointConstraint>(br.a, a0);
  br.problem.add<FixedPointConstraint>(br.b, b0);
  br.problem.add<DistanceConstraint>(br.c, br.a, legLength);
  br.problem.add<DistanceConstraint>(br.c, br.b, legLength);
  br.problem.add<PointOnPlaneConstraint>(
      br.c, Plane3::fromPointNormal(Point3{0, 0, 0}, Vector3{0, 0, 1}));
  return br;
}

bool runCase(const std::string& title, const Point3& initialC, const Point3& expectedC) {
  const Bracket br = buildBracket(initialC);
  const SolverResult result = solve(br.problem);

  std::cout << "=== " << title << " ===\n"
            << std::fixed << std::setprecision(6) << "Initial guess for C : (" << initialC.x
            << ", " << initialC.y << ", " << initialC.z << ")\n";
  std::cout << formatReport(br.problem, result);

  const Point3& c = result.points[br.c.index];
  const double branchError = distance(c, expectedC);
  const double anchorA = distance(result.points[br.a.index], Point3{0, 0, 0});
  const double anchorB = distance(result.points[br.b.index], Point3{100, 0, 0});
  const bool ok = result.converged() && branchError < kCheckTolerance &&
                  anchorA < kCheckTolerance && anchorB < kCheckTolerance &&
                  result.analysis.available && result.analysis.localDegreesOfFreedom == 0;

  std::cout << std::scientific << std::setprecision(3) << "Distance of C to expected solution: "
            << branchError << "\nAnchor drift: A " << anchorA << ", B " << anchorB << "\n"
            << (ok ? "CHECK PASSED" : "CHECK FAILED") << "\n\n";
  return ok;
}

}  // namespace

int main(int argc, char** argv) {
  bool ok = true;
  if (argc == 4) {
    const Point3 guess{std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
    const Point3 expected{50.0, guess.y >= 0.0 ? 50.0 : -50.0, 0.0};
    ok = runCase("custom starting guess", guess, expected);
  } else if (argc == 1) {
    ok = runCase("upper branch (guess y > 0)", Point3{40.0, 30.0, 10.0}, Point3{50, 50, 0}) && ok;
    ok = runCase("lower branch (guess y < 0)", Point3{40.0, -30.0, 10.0}, Point3{50, -50, 0}) && ok;
  } else {
    std::cerr << "usage: bracket_assembly [x y z]\n";
    return 2;
  }
  return ok ? 0 : 1;
}
