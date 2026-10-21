#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "geoconstraint/solver/problem.hpp"

namespace geoconstraint {

/// Residuals of one constraint at the final state.
struct ConstraintResidual {
  std::size_t constraintIndex{0};
  std::string label;
  std::string typeName;
  std::vector<double> values;  ///< `residualCount()` scalar residuals
  double norm{0.0};            ///< Euclidean norm of `values`
};

/// Local, numerical analysis of the constraint Jacobian J (m residuals x n
/// unknowns) at one configuration.
///
/// THIS IS A LOCAL LINEARIZED ESTIMATE, not a proof about the whole nonlinear
/// system: it says how many independent directions of motion are not
/// constrained *to first order* at this configuration. It does not detect
/// multiple discrete solutions, second-order (higher-order) constraints, or
/// global uniqueness, and it depends on `rankTolerance`.
struct SystemAnalysis {
  bool available{false};  ///< false if J could not be computed (see `note`)
  std::string note;       ///< reason when `available` is false
  std::size_t variableCount{0};   ///< n
  std::size_t residualCount{0};   ///< m
  std::size_t rank{0};            ///< numerical rank of J
  std::size_t localDegreesOfFreedom{0};  ///< n - rank
  std::size_t redundantResiduals{0};     ///< m - rank (dependent rows, incl. non-independent components)
  double rankThreshold{0.0};  ///< absolute singular-value cutoff actually used
  double largestSingularValue{0.0};
  double smallestSingularValue{0.0};  ///< over min(m,n) values (0 if m < n)
  std::vector<double> singularValues;  ///< descending, min(m,n) entries
};

/// Rank analysis of `jacobian` by SVD; singular values <= relativeTolerance *
/// sigma_max are treated as zero. Throws InvalidArgumentError for a non-finite
/// matrix or relativeTolerance outside (0,1).
[[nodiscard]] SystemAnalysis analyzeJacobian(const Eigen::MatrixXd& jacobian,
                                             double relativeTolerance);

/// Splits a combined residual vector into per-constraint blocks.
[[nodiscard]] std::vector<ConstraintResidual> perConstraintResiduals(
    const Problem& problem, const Eigen::VectorXd& residuals);

}  // namespace geoconstraint
