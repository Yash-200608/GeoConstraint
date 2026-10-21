#include "geoconstraint/solver/problem.hpp"

#include <string>

namespace geoconstraint {

namespace {
Eigen::Index toIndex(std::size_t n) { return static_cast<Eigen::Index>(n); }
}  // namespace

PointId Problem::addPoint(const Point3& initial, std::string name) {
  const PointId id{initial_.size()};
  if (name.empty()) name = "p" + std::to_string(id.index);
  initial_.push_back(initial);
  pointNames_.push_back(std::move(name));
  return id;
}

std::size_t Problem::addConstraint(std::unique_ptr<Constraint> constraint, std::string label) {
  if (!constraint) throw InvalidArgumentError("Problem::addConstraint: null constraint");
  for (const PointId id : constraint->referencedPoints()) {
    if (id.index >= initial_.size()) {
      throw InvalidArgumentError("Problem::addConstraint: constraint '" +
                                 std::string(constraint->typeName()) + "' references point id " +
                                 std::to_string(id.index) + " but only " +
                                 std::to_string(initial_.size()) + " points exist");
    }
  }
  const std::size_t index = constraints_.size();
  if (label.empty()) label = std::string(constraint->typeName()) + "#" + std::to_string(index);
  offsets_.push_back(residualCount_);
  residualCount_ += constraint->residualCount();
  constraints_.push_back(std::move(constraint));
  constraintLabels_.push_back(std::move(label));
  return index;
}

const std::string& Problem::pointName(PointId id) const {
  if (id.index >= pointNames_.size()) throw InvalidArgumentError("Problem: point id out of range");
  return pointNames_[id.index];
}

const Constraint& Problem::constraint(std::size_t index) const {
  if (index >= constraints_.size()) {
    throw InvalidArgumentError("Problem: constraint index out of range");
  }
  return *constraints_[index];
}

const std::string& Problem::constraintLabel(std::size_t index) const {
  if (index >= constraints_.size()) {
    throw InvalidArgumentError("Problem: constraint index out of range");
  }
  return constraintLabels_[index];
}

std::size_t Problem::residualOffset(std::size_t constraintIndex) const {
  if (constraintIndex >= offsets_.size()) {
    throw InvalidArgumentError("Problem: constraint index out of range");
  }
  return offsets_[constraintIndex];
}

std::vector<std::string> Problem::validate() const {
  std::vector<std::string> issues;
  if (initial_.empty()) issues.emplace_back("problem has no points");
  if (constraints_.empty()) issues.emplace_back("problem has no constraints");
  for (std::size_t i = 0; i < initial_.size(); ++i) {
    if (!initial_[i].isFinite()) {
      issues.push_back("initial position of point '" + pointNames_[i] + "' is not finite");
    }
  }
  return issues;
}

Eigen::VectorXd Problem::encode(std::span<const Point3> points) const {
  if (points.size() != initial_.size()) {
    throw InvalidArgumentError("Problem::encode: wrong number of points");
  }
  Eigen::VectorXd x(toIndex(variableCount()));
  for (std::size_t i = 0; i < points.size(); ++i) {
    x[toIndex(3 * i)] = points[i].x;
    x[toIndex(3 * i + 1)] = points[i].y;
    x[toIndex(3 * i + 2)] = points[i].z;
  }
  return x;
}

std::vector<Point3> Problem::decode(const Eigen::VectorXd& state) const {
  if (state.size() != toIndex(variableCount())) {
    throw InvalidArgumentError("Problem::decode: wrong state size");
  }
  std::vector<Point3> pts(initial_.size());
  for (std::size_t i = 0; i < pts.size(); ++i) {
    pts[i] = {state[toIndex(3 * i)], state[toIndex(3 * i + 1)], state[toIndex(3 * i + 2)]};
  }
  return pts;
}

Eigen::VectorXd Problem::evaluateResiduals(const Eigen::VectorXd& state) const {
  const std::vector<Point3> pts = decode(state);
  Eigen::VectorXd r(toIndex(residualCount_));
  std::vector<double> buffer;
  for (std::size_t c = 0; c < constraints_.size(); ++c) {
    buffer.assign(constraints_[c]->residualCount(), 0.0);
    constraints_[c]->evaluate(pts, buffer);
    for (std::size_t k = 0; k < buffer.size(); ++k) r[toIndex(offsets_[c] + k)] = buffer[k];
  }
  return r;
}

}  // namespace geoconstraint
