#pragma once

#include <span>

#include "geoconstraint/solver/problem.hpp"
#include "geoconstraint/solver/solver_options.hpp"
#include "geoconstraint/solver/solver_result.hpp"

namespace geoconstraint {

/// Minimizes 0.5 ||r(x)||^2 over the point coordinates with a damped
/// least-squares (Levenberg-Marquardt) iteration
///   (J^T J + lambda I) dx = -J^T r,
/// using a finite-difference Jacobian and an LDLT factorization (no explicit
/// inverse). Never throws for bad problems/states: failures are reported through
/// `SolverResult::status` and `message`. The problem is not modified; the
/// starting state is `problem.initialPoints()`.
///
/// Success (`Converged`) is declared only when max|r_i| <= residualTolerance,
/// never merely because the loop ended.
[[nodiscard]] SolverResult solve(const Problem& problem, const SolverOptions& options = {});

/// Finite-difference Jacobian dr/dx at `state` (m x n, row-major residual
/// order as in `Problem`). Step for coordinate j is
/// relativeStep * max(1, |x_j|). Throws if any constraint is invalid at a
/// perturbed state.
[[nodiscard]] Eigen::MatrixXd finiteDifferenceJacobian(
    const Problem& problem, const Eigen::VectorXd& state, double relativeStep = 1e-6,
    JacobianScheme scheme = JacobianScheme::Central);

/// Jacobian-rank / local-DOF analysis at an arbitrary configuration (for
/// example the solution of a previous solve, or a hand-built state). Does not
/// solve anything. Throws GeoConstraintError subclasses if the configuration is
/// invalid for some constraint.
[[nodiscard]] SystemAnalysis analyzeSystem(const Problem& problem, std::span<const Point3> points,
                                           const SolverOptions& options = {});

}  // namespace geoconstraint
