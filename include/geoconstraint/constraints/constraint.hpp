#pragma once

#include <compare>
#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

#include "geoconstraint/math/point3.hpp"

namespace geoconstraint {

/// Stable handle of a point owned by a `Problem` (its index in creation order).
/// Constraints refer to points by id, never by coordinates, so one coordinate
/// is stored exactly once and moves consistently under the solver.
struct PointId {
  std::size_t index{0};
  friend constexpr auto operator<=>(const PointId&, const PointId&) = default;
};

/// A directed line given by two points: the direction vector is
/// `P[to] - P[from]`. Only the *direction* matters to the constraints that take
/// a `DirectionRef`; the position of the line is irrelevant.
struct DirectionRef {
  PointId from;
  PointId to;
};

/// Interface of a scalar-residual geometric constraint.
///
/// A constraint is satisfied when all of its residuals are zero. Residual
/// definitions are given in each concrete class and in
/// docs/mathematical_foundations.md. `residualCount()` is a fixed property of
/// the constraint type, independent of the geometry, so the combined residual
/// vector layout is deterministic.
class Constraint {
 public:
  Constraint() = default;
  Constraint(const Constraint&) = default;
  Constraint& operator=(const Constraint&) = default;
  Constraint(Constraint&&) = default;
  Constraint& operator=(Constraint&&) = default;
  virtual ~Constraint() = default;

  /// Stable human-readable type name, e.g. "Distance".
  [[nodiscard]] virtual std::string_view typeName() const noexcept = 0;

  /// Number of scalar residuals this constraint produces.
  [[nodiscard]] virtual std::size_t residualCount() const noexcept = 0;

  /// All point ids this constraint reads (may contain duplicates).
  [[nodiscard]] virtual std::vector<PointId> referencedPoints() const = 0;

  /// Writes `residualCount()` residuals into `residuals`.
  /// `points` is indexed by `PointId::index`.
  /// Throws InvalidArgumentError if an id is out of range or the output span has
  /// the wrong size, and DegenerateGeometryError if a direction used by the
  /// constraint is zero-length / non-finite (such a direction must never
  /// satisfy a parallelism/perpendicularity/angle constraint by accident).
  virtual void evaluate(std::span<const Point3> points, std::span<double> residuals) const = 0;

  /// Convenience wrapper allocating the output.
  [[nodiscard]] std::vector<double> residuals(std::span<const Point3> points) const;

 protected:
  /// Bounds-checked point lookup.
  [[nodiscard]] static const Point3& pointAt(std::span<const Point3> points, PointId id);
  /// Unit direction of a `DirectionRef`; throws DegenerateGeometryError if its
  /// length is not above `Tolerances::zeroLength` or a coordinate is non-finite.
  [[nodiscard]] static Vector3 unitDirection(std::span<const Point3> points, DirectionRef ref);
  /// Throws InvalidArgumentError unless `residuals.size() == expected`.
  static void requireSize(std::span<double> residuals, std::size_t expected);
};

}  // namespace geoconstraint
