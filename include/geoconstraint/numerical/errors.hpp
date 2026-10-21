#pragma once

#include <stdexcept>

namespace geoconstraint {

/// Base class of every exception thrown by GeoConstraint.
class GeoConstraintError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

/// An argument is outside its documented domain (non-finite number, negative
/// distance, index out of range, matrix that is not a rotation, ...).
class InvalidArgumentError : public GeoConstraintError {
 public:
  using GeoConstraintError::GeoConstraintError;
};

/// The input is numerically well-formed but degenerate for the requested
/// operation (zero-length direction, collinear points defining a plane, ...).
class DegenerateGeometryError : public GeoConstraintError {
 public:
  using GeoConstraintError::GeoConstraintError;
};

}  // namespace geoconstraint
