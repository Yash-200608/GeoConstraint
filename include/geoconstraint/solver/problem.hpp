#pragma once

#include <cstddef>
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "geoconstraint/constraints/constraint.hpp"

namespace geoconstraint {

/// A constraint system: points (the unknowns, each with an initial position)
/// and constraints over them.
///
/// State layout: x = [x0 y0 z0 x1 y1 z1 ...], i.e. point k occupies entries
/// 3k..3k+2. The combined residual vector concatenates the constraints'
/// residuals in the order they were added; `residualOffset(i)` gives the start
/// of constraint i. Both layouts are deterministic.
///
/// Every point is an unknown. To "fix" a point, add a `FixedPointConstraint`.
class Problem {
 public:
  /// Adds a point with its initial position (the solver's starting guess).
  /// Non-finite coordinates are accepted here and reported by `validate()` /
  /// the solver as InvalidInput, so a bad initial state yields a structured
  /// result rather than an exception.
  PointId addPoint(const Point3& initial, std::string name = {});

  /// Takes ownership of `constraint`. Throws InvalidArgumentError if it is null
  /// or references a point id that does not exist yet. Returns its index.
  std::size_t addConstraint(std::unique_ptr<Constraint> constraint, std::string label = {});

  /// Constructs a constraint of type C in place.
  template <class C, class... Args>
  std::size_t add(Args&&... args) {
    return addConstraint(std::make_unique<C>(std::forward<Args>(args)...));
  }

  [[nodiscard]] std::size_t pointCount() const noexcept { return initial_.size(); }
  [[nodiscard]] std::size_t constraintCount() const noexcept { return constraints_.size(); }
  [[nodiscard]] std::size_t variableCount() const noexcept { return 3 * initial_.size(); }
  [[nodiscard]] std::size_t residualCount() const noexcept { return residualCount_; }

  [[nodiscard]] std::span<const Point3> initialPoints() const noexcept { return initial_; }
  [[nodiscard]] const std::string& pointName(PointId id) const;
  [[nodiscard]] const Constraint& constraint(std::size_t index) const;
  [[nodiscard]] const std::string& constraintLabel(std::size_t index) const;
  [[nodiscard]] std::size_t residualOffset(std::size_t constraintIndex) const;

  /// Structural/initial-state problems that make solving meaningless (no
  /// points, no constraints, non-finite initial coordinates). Empty if valid.
  [[nodiscard]] std::vector<std::string> validate() const;

  /// Flattens points into the state vector. Throws if the span has the wrong size.
  [[nodiscard]] Eigen::VectorXd encode(std::span<const Point3> points) const;
  /// Inverse of `encode`. Throws InvalidArgumentError on a wrong-sized state.
  [[nodiscard]] std::vector<Point3> decode(const Eigen::VectorXd& state) const;

  /// Combined residual vector r(x). Throws GeoConstraintError subclasses if any
  /// constraint is invalid at x (e.g. a zero-length direction).
  [[nodiscard]] Eigen::VectorXd evaluateResiduals(const Eigen::VectorXd& state) const;

 private:
  std::vector<Point3> initial_;
  std::vector<std::string> pointNames_;
  std::vector<std::unique_ptr<Constraint>> constraints_;
  std::vector<std::string> constraintLabels_;
  std::vector<std::size_t> offsets_;
  std::size_t residualCount_{0};
};

}  // namespace geoconstraint
