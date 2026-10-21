#include "geoconstraint/solver/solver.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace geoconstraint {

namespace {

using Eigen::Index;
using Eigen::MatrixXd;
using Eigen::VectorXd;

constexpr auto kInf = Eigen::Infinity;

std::size_t evaluationsPerJacobian(const Problem& problem, JacobianScheme scheme) {
  return problem.variableCount() * (scheme == JacobianScheme::Central ? 2 : 1);
}

}  // namespace

std::optional<std::string> SolverOptions::validate() const {
  const auto positive = [](double v) { return std::isfinite(v) && v > 0.0; };
  if (maxIterations < 1) return "maxIterations must be >= 1";
  if (!positive(residualTolerance)) return "residualTolerance must be finite and > 0";
  if (!(std::isfinite(stepAbsTolerance) && stepAbsTolerance >= 0.0) ||
      !(std::isfinite(stepRelTolerance) && stepRelTolerance >= 0.0)) {
    return "step tolerances must be finite and >= 0";
  }
  if (!(std::isfinite(gradientTolerance) && gradientTolerance >= 0.0)) {
    return "gradientTolerance must be finite and >= 0";
  }
  if (!positive(initialDamping)) return "initialDamping must be finite and > 0";
  if (!(std::isfinite(dampingDecrease) && dampingDecrease > 0.0 && dampingDecrease < 1.0)) {
    return "dampingDecrease must be in (0,1)";
  }
  if (!(std::isfinite(dampingIncrease) && dampingIncrease > 1.0)) {
    return "dampingIncrease must be > 1";
  }
  if (!positive(minDamping) || !positive(maxDamping) || minDamping > maxDamping) {
    return "damping bounds must satisfy 0 < minDamping <= maxDamping";
  }
  if (initialDamping < minDamping || initialDamping > maxDamping) {
    return "initialDamping must lie within [minDamping, maxDamping]";
  }
  if (!(std::isfinite(minGainRatio) && minGainRatio >= 0.0 && minGainRatio < 1.0)) {
    return "minGainRatio must be in [0,1)";
  }
  if (!positive(finiteDifferenceStep)) return "finiteDifferenceStep must be finite and > 0";
  if (!(std::isfinite(rankTolerance) && rankTolerance > 0.0 && rankTolerance < 1.0)) {
    return "rankTolerance must be in (0,1)";
  }
  return std::nullopt;
}

MatrixXd finiteDifferenceJacobian(const Problem& problem, const VectorXd& state,
                                  double relativeStep, JacobianScheme scheme) {
  if (!(std::isfinite(relativeStep) && relativeStep > 0.0)) {
    throw InvalidArgumentError("finiteDifferenceJacobian: step must be finite and > 0");
  }
  const Index n = state.size();
  const VectorXd r0 = scheme == JacobianScheme::Forward ? problem.evaluateResiduals(state)
                                                        : VectorXd();
  MatrixXd jac(static_cast<Index>(problem.residualCount()), n);
  VectorXd probe = state;
  for (Index j = 0; j < n; ++j) {
    // Scale the step with the coordinate magnitude, then use the step that was
    // actually representable in floating point (x+h)-x to avoid rounding bias.
    const double nominal = relativeStep * std::max(1.0, std::abs(state[j]));
    probe[j] = state[j] + nominal;
    const double hPlus = probe[j] - state[j];
    const VectorXd rPlus = problem.evaluateResiduals(probe);
    if (scheme == JacobianScheme::Central) {
      probe[j] = state[j] - nominal;
      const double hMinus = state[j] - probe[j];
      const VectorXd rMinus = problem.evaluateResiduals(probe);
      jac.col(j) = (rPlus - rMinus) / (hPlus + hMinus);
    } else {
      jac.col(j) = (rPlus - r0) / hPlus;
    }
    probe[j] = state[j];
  }
  return jac;
}

SystemAnalysis analyzeSystem(const Problem& problem, std::span<const Point3> points,
                             const SolverOptions& options) {
  if (auto err = options.validate()) throw InvalidArgumentError("analyzeSystem: " + *err);
  const VectorXd x = problem.encode(points);
  const MatrixXd jac = finiteDifferenceJacobian(problem, x, options.finiteDifferenceStep,
                                                options.jacobianScheme);
  return analyzeJacobian(jac, options.rankTolerance);
}

SolverResult solve(const Problem& problem, const SolverOptions& opt) {
  SolverResult res;
  res.points.assign(problem.initialPoints().begin(), problem.initialPoints().end());
  res.finalDamping = opt.initialDamping;

  const auto fail = [&res](SolverStatus status, std::string message) {
    res.status = status;
    res.message = std::move(message);
    res.analysis.note = "not computed: the solve did not start";
    return res;
  };

  if (const auto err = opt.validate()) return fail(SolverStatus::InvalidInput, "invalid options: " + *err);
  if (const auto issues = problem.validate(); !issues.empty()) {
    std::string msg = "invalid problem: ";
    for (std::size_t i = 0; i < issues.size(); ++i) msg += (i ? "; " : "") + issues[i];
    return fail(SolverStatus::InvalidInput, std::move(msg));
  }

  VectorXd x = problem.encode(problem.initialPoints());
  VectorXd r;
  try {
    r = problem.evaluateResiduals(x);
  } catch (const GeoConstraintError& e) {
    return fail(SolverStatus::InvalidInput,
                std::string("constraints cannot be evaluated at the initial state: ") + e.what());
  }
  if (!r.allFinite()) {
    return fail(SolverStatus::InvalidInput, "non-finite residual at the initial state");
  }
  res.residualEvaluations = 1;
  res.initialResidualNorm = r.norm();
  res.costHistory.push_back(0.5 * r.squaredNorm());

  // Records the final state, residual breakdown and Jacobian rank/DOF estimate.
  const auto finish = [&](SolverStatus status, std::string message) {
    res.status = status;
    res.message = std::move(message);
    res.points = problem.decode(x);
    res.finalResidualNorm = r.norm();
    res.finalResidualMax = r.lpNorm<kInf>();
    res.finalCost = 0.5 * r.squaredNorm();
    res.constraintResiduals = perConstraintResiduals(problem, r);
    try {
      const MatrixXd jac = finiteDifferenceJacobian(problem, x, opt.finiteDifferenceStep,
                                                    opt.jacobianScheme);
      res.residualEvaluations += evaluationsPerJacobian(problem, opt.jacobianScheme);
      res.analysis = analyzeJacobian(jac, opt.rankTolerance);
      res.finalGradientNorm = (jac.transpose() * r).lpNorm<kInf>();
    } catch (const GeoConstraintError& e) {
      res.analysis = SystemAnalysis{};
      res.analysis.note = e.what();
    }
    return res;
  };

  double lambda = opt.initialDamping;

  for (;;) {
    res.finalDamping = lambda;
    if (r.lpNorm<kInf>() <= opt.residualTolerance) {
      return finish(SolverStatus::Converged, "all residuals within residualTolerance");
    }
    if (res.iterations >= static_cast<std::size_t>(opt.maxIterations)) {
      return finish(SolverStatus::MaxIterations,
                    "iteration limit reached before the residual tolerance was met");
    }

    // Steps 3-4: finite-difference Jacobian, gradient and Gauss-Newton Hessian.
    MatrixXd jac;
    try {
      jac = finiteDifferenceJacobian(problem, x, opt.finiteDifferenceStep, opt.jacobianScheme);
    } catch (const GeoConstraintError& e) {
      return finish(SolverStatus::NumericalFailure,
                    std::string("Jacobian evaluation failed: ") + e.what());
    }
    res.residualEvaluations += evaluationsPerJacobian(problem, opt.jacobianScheme);
    if (!jac.allFinite()) {
      return finish(SolverStatus::NumericalFailure, "non-finite entries in the Jacobian");
    }
    const VectorXd grad = jac.transpose() * r;
    if (grad.lpNorm<kInf>() <= opt.gradientTolerance) {
      return finish(SolverStatus::Stagnated,
                    "stationary point of the least-squares objective with non-zero residual "
                    "(inconsistent constraints or a local minimum)");
    }
    const MatrixXd jtj = jac.transpose() * jac;
    const double cost = 0.5 * r.squaredNorm();

    // Steps 5-8: damped solve; raise damping until a step is accepted.
    bool accepted = false;
    bool lastFailureWasLinearSolve = false;
    while (!accepted) {
      if (lambda > opt.maxDamping) {
        return finish(lastFailureWasLinearSolve ? SolverStatus::NumericalFailure
                                                : SolverStatus::Stagnated,
                      lastFailureWasLinearSolve
                          ? "linear solve failed for every damping value"
                          : "no damped step reduces the objective (damping exceeded maximum); "
                            "stationary point, inconsistent constraints, or Jacobian noise");
      }
      MatrixXd lhs = jtj;
      lhs.diagonal().array() += lambda;
      const Eigen::LDLT<MatrixXd> ldlt(lhs);
      VectorXd step;
      bool solved = ldlt.info() == Eigen::Success;
      if (solved) {
        step = ldlt.solve(-grad);
        solved = step.allFinite();
      }
      if (!solved) {
        lastFailureWasLinearSolve = true;
        lambda *= opt.dampingIncrease;
        ++res.rejectedSteps;
        continue;
      }
      lastFailureWasLinearSolve = false;

      const VectorXd trial = x + step;
      VectorXd trialResidual;
      bool valid = true;
      try {
        trialResidual = problem.evaluateResiduals(trial);
        ++res.residualEvaluations;
        valid = trialResidual.allFinite();
      } catch (const GeoConstraintError&) {
        valid = false;  // trial step entered a degenerate region: reject it
      }

      if (valid) {
        // Gain ratio: actual reduction vs. reduction predicted by the damped
        // quadratic model, 0.5 * step^T (lambda*step - grad).
        const double actual = cost - 0.5 * trialResidual.squaredNorm();
        const double predicted = 0.5 * step.dot(lambda * step - grad);
        accepted = actual > 0.0 && predicted > 0.0 && actual / predicted >= opt.minGainRatio;
      }

      if (accepted) {
        const double stepNorm = step.norm();
        x = trial;
        r = trialResidual;
        lambda = std::max(lambda * opt.dampingDecrease, opt.minDamping);
        ++res.iterations;
        res.costHistory.push_back(0.5 * r.squaredNorm());
        res.finalDamping = lambda;
        const bool stepTiny =
            stepNorm <= opt.stepAbsTolerance + opt.stepRelTolerance * x.norm();
        if (stepTiny && r.lpNorm<kInf>() > opt.residualTolerance) {
          return finish(SolverStatus::Stagnated,
                        "step size below tolerance while the residual is still too large "
                        "(inconsistent constraints or a local minimum)");
        }
      } else {
        lambda *= opt.dampingIncrease;
        ++res.rejectedSteps;
      }
    }
  }
}

}  // namespace geoconstraint
