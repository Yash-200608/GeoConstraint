// Shows how the solver reports systems that are NOT a clean unique solution:
// underconstrained, inconsistent, iteration-limited and invalid problems.
// The program exits non-zero if a scenario does not report the expected status.

#include <iostream>
#include <limits>
#include <string>

#include "geoconstraint/constraints/coincidence_constraint.hpp"
#include "geoconstraint/constraints/distance_constraint.hpp"
#include "geoconstraint/constraints/fixed_point_constraint.hpp"
#include "geoconstraint/constraints/parallel_constraint.hpp"
#include "geoconstraint/solver/solver.hpp"

using namespace geoconstraint;

namespace {

int failures = 0;

void expect(const std::string& title, const Problem& problem, const SolverResult& result,
            SolverStatus wanted) {
  std::cout << "=== " << title << " ===\n" << formatReport(problem, result);
  const bool ok = result.status == wanted;
  std::cout << (ok ? "EXPECTED: " : "UNEXPECTED, wanted: ") << toString(wanted) << "\n\n";
  if (!ok) ++failures;
}

}  // namespace

int main() {
  {  // Converges, yet two directions of motion remain free.
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}, "A");
    const PointId c = p.addPoint({1, 2, 2}, "C");
    p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
    p.add<DistanceConstraint>(c, a, 5.0);
    const SolverResult r = solve(p);
    expect("Underdetermined: C on a sphere around A (expect local DOF = 2)", p, r,
           SolverStatus::Converged);
    if (r.analysis.localDegreesOfFreedom != 2) {
      std::cout << "UNEXPECTED local DOF " << r.analysis.localDegreesOfFreedom << "\n";
      ++failures;
    }
  }
  {  // Contradictory requirements: no point is 1 unit from both A and B, 10 apart.
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}, "A");
    const PointId b = p.addPoint({10, 0, 0}, "B");
    const PointId c = p.addPoint({5, 1, 0}, "C");
    p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
    p.add<FixedPointConstraint>(b, Point3{10, 0, 0});
    p.add<DistanceConstraint>(c, a, 1.0);
    p.add<DistanceConstraint>(c, b, 1.0);
    expect("Inconsistent constraints (no solution exists)", p, solve(p), SolverStatus::Stagnated);
  }
  {  // Iteration cap.
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}, "A");
    const PointId c = p.addPoint({1, 1, 1}, "C");
    p.add<FixedPointConstraint>(a, Point3{0, 0, 0});
    p.add<DistanceConstraint>(c, a, 1000.0);
    SolverOptions options;
    options.maxIterations = 1;
    expect("Iteration limit of 1", p, solve(p, options), SolverStatus::MaxIterations);
  }
  {  // Invalid: a zero-length direction in a parallel constraint.
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}, "A");
    const PointId b = p.addPoint({0, 0, 0}, "B");  // coincides with A -> no direction
    const PointId c = p.addPoint({0, 0, 0}, "C");
    const PointId d = p.addPoint({1, 0, 0}, "D");
    p.add<ParallelConstraint>(DirectionRef{a, b}, DirectionRef{c, d});
    expect("Invalid input: zero-length direction", p, solve(p), SolverStatus::InvalidInput);
  }
  {  // Invalid: non-finite coordinate.
    Problem p;
    const PointId a = p.addPoint({0, 0, 0}, "A");
    const PointId b = p.addPoint({std::numeric_limits<double>::quiet_NaN(), 0, 0}, "B");
    p.add<CoincidenceConstraint>(a, b);
    expect("Invalid input: NaN coordinate", p, solve(p), SolverStatus::InvalidInput);
  }

  std::cout << (failures == 0 ? "All scenarios reported as expected\n"
                              : "Some scenarios did NOT report as expected\n");
  return failures == 0 ? 0 : 1;
}
