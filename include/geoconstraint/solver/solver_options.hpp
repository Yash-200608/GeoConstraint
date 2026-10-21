#pragma once

#include <optional>
#include <string>

namespace geoconstraint {

/// Finite-difference scheme for the Jacobian.
enum class JacobianScheme {
  Forward,  ///< (r(x+h) - r(x)) / h            : O(h) truncation, 1 eval/column
  Central   ///< (r(x+h) - r(x-h)) / (2h)       : O(h^2) truncation, 2 evals/column
};

/// Configuration of the Levenberg-Marquardt solver. The defaults suit
/// coordinates of order 1e-2 ... 1e4 model units; see docs/solver_design.md for
/// how each option enters the algorithm and its limitations.
struct SolverOptions {
  /// Hard cap on *accepted* iterations.
  int maxIterations = 100;

  /// Success criterion: max_i |r_i| <= residualTolerance. ABSOLUTE, in the units
  /// of the residuals (lengths for distance/point constraints, dimensionless for
  /// direction constraints). Chosen so the final geometry is accurate to
  /// ~1e-9 model units for unit-scale problems.
  double residualTolerance = 1e-9;

  /// Stagnation: an accepted step with ||dx|| <= stepAbsTolerance +
  /// stepRelTolerance * ||x|| while the residual criterion is not met.
  double stepAbsTolerance = 1e-12;
  double stepRelTolerance = 1e-12;

  /// Stagnation: ||J^T r||_inf <= gradientTolerance while the residual
  /// criterion is not met (stationary point of the least-squares objective).
  double gradientTolerance = 1e-10;

  /// Initial Levenberg-Marquardt damping lambda (adds lambda*I to J^T J).
  double initialDamping = 1e-3;
  /// Damping multiplier after an accepted step (in (0,1)) / a rejected step (>1).
  double dampingDecrease = 0.3;
  double dampingIncrease = 4.0;
  /// Damping is clamped below by minDamping; exceeding maxDamping on a rejected
  /// step terminates the solve (no damped step improves the objective).
  double minDamping = 1e-14;
  double maxDamping = 1e12;

  /// A step is accepted only if the gain ratio
  ///   rho = (actual cost reduction) / (reduction predicted by the model)
  /// is at least this value (and the actual reduction is positive).
  double minGainRatio = 1e-4;

  /// Finite-difference step for coordinate j: h_j = step * max(1, |x_j|).
  /// Central differences with 1e-6 balance truncation (~h^2) and rounding
  /// (~eps/h) errors to roughly 1e-10 relative accuracy.
  double finiteDifferenceStep = 1e-6;
  JacobianScheme jacobianScheme = JacobianScheme::Central;

  /// Numerical-rank threshold for the diagnostics: singular values
  /// <= rankTolerance * sigma_max count as zero. Must be well above the
  /// finite-difference noise (~1e-10); see docs/numerical_robustness.md.
  double rankTolerance = 1e-8;

  /// Returns a description of the first invalid option, or std::nullopt.
  [[nodiscard]] std::optional<std::string> validate() const;
};

}  // namespace geoconstraint
