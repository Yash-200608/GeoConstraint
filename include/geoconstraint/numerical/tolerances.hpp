#pragma once

#include <algorithm>
#include <cmath>

namespace geoconstraint {

/// Central numerical policy. Every tolerance used by the geometry layer lives
/// here; there is no other hidden threshold in the geometry code.
///
/// Intended coordinate scale: model units of roughly 1e-3 ... 1e6 (e.g. mm for
/// mechanical parts). Absolute tolerances assume that range; the `relative`
/// term widens length comparisons for large coordinates, because doubles near
/// 1e6 only resolve ~1e-10 and rounding grows with magnitude.
struct Tolerances {
  /// A vector shorter than this (absolute, model units) is "zero" and cannot
  /// define a direction. Far below any meaningful geometry, far above denormals.
  double zeroLength = 1e-12;
  /// Absolute length tolerance for coincidence / containment tests.
  double length = 1e-9;
  /// Relative length tolerance, multiplied by the magnitude (scale) of the
  /// operands. Effective length tolerance = length + relative * scale.
  double relative = 1e-9;
  /// Dimensionless angular tolerance: sin(angle) for parallelism, cos(angle)
  /// for perpendicularity (both evaluated on normalized directions). 1e-9 rad
  /// is ~6e-8 degrees.
  double angular = 1e-9;
};

/// The default policy. `constexpr`, so there is no mutable global state; pass a
/// different `Tolerances` explicitly to override.
inline constexpr Tolerances kDefaultTolerances{};

/// Effective length tolerance for operands of magnitude `scale`.
[[nodiscard]] inline double lengthTolerance(const Tolerances& tol, double scale) noexcept {
  return tol.length + tol.relative * std::abs(scale);
}

/// |a-b| <= length + relative * max(|a|,|b|). Non-finite inputs never compare equal.
[[nodiscard]] inline bool approxEqual(double a, double b,
                                      const Tolerances& tol = kDefaultTolerances) noexcept {
  if (!std::isfinite(a) || !std::isfinite(b)) return false;
  return std::abs(a - b) <= lengthTolerance(tol, std::max(std::abs(a), std::abs(b)));
}

}  // namespace geoconstraint
