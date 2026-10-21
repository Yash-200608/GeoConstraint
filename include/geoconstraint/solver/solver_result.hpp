#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "geoconstraint/solver/diagnostics.hpp"

namespace geoconstraint {

/// How a solve ended. Only `Converged` means the constraints are satisfied.
enum class SolverStatus {
  Converged,         ///< max |r_i| <= residualTolerance
  MaxIterations,     ///< iteration limit reached before convergence or stagnation
  Stagnated,         ///< no further progress possible but residual still too large:
                     ///< inconsistent constraints or a local minimum of the objective
  InvalidInput,      ///< bad problem/options/initial state, or a constraint is invalid
                     ///< (e.g. zero-length direction) at the starting configuration
  NumericalFailure   ///< non-finite values, failed linear solves, or the Jacobian
                     ///< could not be evaluated mid-solve
};

[[nodiscard]] std::string_view toString(SolverStatus status) noexcept;

/// Everything the solver reports. `points` is always filled: for Converged it
/// is the solution; for other statuses it is the best (lowest-cost) state
/// reached, which is *not* a valid solution.
struct SolverResult {
  SolverStatus status{SolverStatus::InvalidInput};
  std::string message;  ///< human-readable reason for `status`

  std::vector<Point3> points;  ///< final state, indexed by PointId

  std::size_t iterations{0};       ///< accepted LM iterations
  std::size_t rejectedSteps{0};    ///< trial steps rejected (incl. failed linear solves)
  std::size_t residualEvaluations{0};  ///< total calls to r(x), incl. finite differences

  double initialResidualNorm{0.0};  ///< ||r||_2 at the starting state
  double finalResidualNorm{0.0};    ///< ||r||_2 at `points`
  double finalResidualMax{0.0};     ///< ||r||_inf at `points`
  double finalCost{0.0};            ///< 0.5 ||r||^2
  double finalGradientNorm{0.0};    ///< ||J^T r||_inf at `points` (0 if unavailable)
  double finalDamping{0.0};         ///< LM lambda at termination

  std::vector<double> costHistory;  ///< 0.5||r||^2 at the start and after each accepted step
  std::vector<ConstraintResidual> constraintResiduals;  ///< at `points`
  SystemAnalysis analysis;  ///< local Jacobian rank / DOF estimate at `points`

  [[nodiscard]] bool converged() const noexcept { return status == SolverStatus::Converged; }
};

/// Multi-line, human-readable report (status, residuals per constraint, point
/// coordinates with the problem's point names, Jacobian rank and local DOF).
[[nodiscard]] std::string formatReport(const Problem& problem, const SolverResult& result);

}  // namespace geoconstraint
