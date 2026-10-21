#include "geoconstraint/solver/diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "geoconstraint/solver/solver_result.hpp"

namespace geoconstraint {

namespace {
Eigen::Index toIndex(std::size_t n) { return static_cast<Eigen::Index>(n); }
}  // namespace

SystemAnalysis analyzeJacobian(const Eigen::MatrixXd& jacobian, double relativeTolerance) {
  if (!(relativeTolerance > 0.0 && relativeTolerance < 1.0)) {
    throw InvalidArgumentError("analyzeJacobian: relativeTolerance must be in (0,1)");
  }
  if (!jacobian.allFinite()) throw InvalidArgumentError("analyzeJacobian: non-finite Jacobian");

  SystemAnalysis a;
  a.variableCount = static_cast<std::size_t>(jacobian.cols());
  a.residualCount = static_cast<std::size_t>(jacobian.rows());
  const Eigen::BDCSVD<Eigen::MatrixXd> svd(jacobian);  // singular values only
  const Eigen::VectorXd sv = svd.singularValues();      // descending
  a.singularValues.assign(sv.data(), sv.data() + sv.size());
  a.largestSingularValue = sv.size() > 0 ? sv[0] : 0.0;
  a.smallestSingularValue = sv.size() > 0 ? sv[sv.size() - 1] : 0.0;
  a.rankThreshold = relativeTolerance * a.largestSingularValue;
  a.rank = static_cast<std::size_t>(
      std::count_if(a.singularValues.begin(), a.singularValues.end(),
                    [&](double s) { return s > a.rankThreshold; }));
  a.localDegreesOfFreedom = a.variableCount - a.rank;
  a.redundantResiduals = a.residualCount - a.rank;
  a.available = true;
  return a;
}

std::vector<ConstraintResidual> perConstraintResiduals(const Problem& problem,
                                                       const Eigen::VectorXd& residuals) {
  if (residuals.size() != toIndex(problem.residualCount())) {
    throw InvalidArgumentError("perConstraintResiduals: wrong residual vector size");
  }
  std::vector<ConstraintResidual> out;
  out.reserve(problem.constraintCount());
  for (std::size_t c = 0; c < problem.constraintCount(); ++c) {
    ConstraintResidual cr;
    cr.constraintIndex = c;
    cr.label = problem.constraintLabel(c);
    cr.typeName = std::string(problem.constraint(c).typeName());
    const std::size_t off = problem.residualOffset(c);
    const std::size_t count = problem.constraint(c).residualCount();
    double sumSq = 0.0;
    for (std::size_t k = 0; k < count; ++k) {
      const double v = residuals[toIndex(off + k)];
      cr.values.push_back(v);
      sumSq += v * v;
    }
    cr.norm = std::sqrt(sumSq);
    out.push_back(std::move(cr));
  }
  return out;
}

std::string_view toString(SolverStatus status) noexcept {
  switch (status) {
    case SolverStatus::Converged: return "Converged";
    case SolverStatus::MaxIterations: return "MaxIterations";
    case SolverStatus::Stagnated: return "Stagnated";
    case SolverStatus::InvalidInput: return "InvalidInput";
    case SolverStatus::NumericalFailure: return "NumericalFailure";
  }
  return "Unknown";
}

std::string formatReport(const Problem& problem, const SolverResult& r) {
  std::ostringstream os;
  os << std::setprecision(12);
  os << "Status        : " << toString(r.status) << " - " << r.message << "\n";
  os << "Iterations    : " << r.iterations << " accepted, " << r.rejectedSteps
     << " rejected steps, " << r.residualEvaluations << " residual evaluations\n";
  os << std::scientific << std::setprecision(3);
  os << "Residual      : |r|_2 = " << r.finalResidualNorm << " (initial " << r.initialResidualNorm
     << "), |r|_inf = " << r.finalResidualMax << "\n";
  os << "Gradient      : |J^T r|_inf = " << r.finalGradientNorm << ", final damping = "
     << r.finalDamping << "\n";
  os << std::fixed << std::setprecision(9);
  os << "Points:\n";
  for (std::size_t i = 0; i < r.points.size() && i < problem.pointCount(); ++i) {
    const Point3& p = r.points[i];
    os << "  " << problem.pointName(PointId{i}) << " = (" << p.x << ", " << p.y << ", " << p.z
       << ")\n";
  }
  os << "Constraint residuals:\n" << std::scientific << std::setprecision(3);
  for (const ConstraintResidual& c : r.constraintResiduals) {
    os << "  [" << c.constraintIndex << "] " << c.label << " (" << c.typeName << "): r = [";
    for (std::size_t k = 0; k < c.values.size(); ++k) os << (k ? ", " : "") << c.values[k];
    os << "], |r| = " << c.norm << "\n";
  }
  const SystemAnalysis& a = r.analysis;
  if (a.available) {
    os << "Jacobian (local numerical estimate): " << a.residualCount << " x " << a.variableCount
       << ", rank = " << a.rank << " (threshold " << a.rankThreshold << "), local DOF = "
       << a.localDegreesOfFreedom << ", redundant residual rows = " << a.redundantResiduals
       << ", sigma_max = " << a.largestSingularValue << ", sigma_min = "
       << a.smallestSingularValue << "\n";
  } else {
    os << "Jacobian analysis unavailable: " << a.note << "\n";
  }
  return os.str();
}

}  // namespace geoconstraint
