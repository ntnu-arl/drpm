#include <boost/math/distributions/normal.hpp>
#include <cmath>
#include "Eigen/Eigenvalues"

namespace degeneracy {

template <typename T>
using VectorVector3 = std::vector<Eigen::Matrix<T, 3, 1>>;

template <typename T>
using VectorMatrix3 = std::vector<Eigen::Matrix<T, 3, 3>>;

template <typename T>
inline Eigen::Matrix<T, 3, 3> VectorToSkew(const Eigen::Matrix<T, 3, 1>& vector) {
  Eigen::Matrix<T, 3, 3> skew;
  skew << 0, -vector.z(), vector.y(), vector.z(), 0, -vector.x(), -vector.y(), vector.x(), 0;
  return skew;
}

template <typename T, typename Q>
auto ComputeNoiseEstimate(const VectorVector3<Q>& points, const VectorVector3<Q>& normals, const std::vector<Q>& weights, const VectorMatrix3<T>& normal_covariances, const Eigen::Matrix<T, 6, 6>& U, const T& stdevPoints) {
  using Vector3 = Eigen::Matrix<T, 3, 1>;
  using Matrix3 = Eigen::Matrix<T, 3, 3>;
  using Vector6 = Eigen::Matrix<T, 6, 1>;
  using Matrix6 = Eigen::Matrix<T, 6, 6>;

  Matrix6 mean = Matrix6::Zero();
  Vector6 variance = Vector6::Zero();

  const size_t nPoints = points.size();

  for (size_t i = 0; i < nPoints; i++) {
    const Vector3 point = points[i].template cast<T>();
    const Vector3 normal = normals[i].template cast<T>();
    const Matrix3 nx = VectorToSkew<T>(normal);
    const Matrix3 px = VectorToSkew<T>(point);
    const T w = weights[i];

    // Coefficient matrix for epsilon and eta
    Matrix6 B = Matrix6::Zero();
    B.block(0, 0, 3, 3) = -nx;
    B.block(0, 3, 3, 3) = px * nx;
    B.block(3, 3, 3, 3) = nx;

    // Covariance matrix for epsilon and eta
    Matrix6 N = Matrix6::Zero();
    N.block(0, 0, 3, 3) = Matrix3::Identity() * std::pow(stdevPoints, 2);
    N.block(3, 3, 3, 3) = normal_covariances[i];

    Matrix6 contribution_to_mean = (B * N * B.transpose()) * w;

    mean.noalias() += contribution_to_mean.eval();

    // v hat weighted by w
    Vector6 v = Vector6::Zero();
    v.head(3) = std::sqrt(w) * px * normal;
    v.tail(3) = std::sqrt(w) * normal;

    // Compute variance in the directions given by U
    for (size_t k = 0; k < 6; k++) {
      const Vector6 u = U.col(k);
      const T a = (u.transpose() * contribution_to_mean * u).value();
      const T b = (u.transpose() * v).value();
      const T contribution_to_variance = 2 * std::pow(a, 2) + 4 * a * std::pow(b, 2);
      variance[k] += contribution_to_variance;
    }
  }

  return std::make_tuple(mean, variance);
}

template <typename T>
Eigen::Matrix<T, 6, 1> ComputeSignalToNoiseProbabilities(const Eigen::Matrix<T, 6, 6>& measured_information_matrix,
                                                         const Eigen::Matrix<T, 6, 6>& estimated_noise_mean,
                                                         const Eigen::Matrix<T, 6, 1>& estimated_noise_variances,
                                                         const Eigen::Matrix<T, 6, 6>& U,
                                                         const T& snr_factor) {
  typedef Eigen::Matrix<T, 6, 1> Vector6;

  Vector6 probabilities = Vector6::Zero();

  for (size_t k = 0; k < 6; k++) {
    const Vector6 u = U.col(k);
    const T measurement = (u.transpose() * measured_information_matrix * u).value();
    const T expected_noise = (u.transpose() * estimated_noise_mean * u).value();
    const T stdev = std::sqrt(estimated_noise_variances[k]);
    const T test_point = measurement / (T(1.0) + snr_factor);

    const bool any_nan = std::isnan(expected_noise) || std::isnan(stdev) || std::isnan(test_point);

    if (!any_nan) {
      const T probability = boost::math::cdf(
          boost::math::normal_distribution<T>(expected_noise, stdev), test_point);

      probabilities[k] = probability;
    } else {
      std::cout << "NaN value in probability calculation - stDev: " << stdev << " test point: " << test_point << " expected noise: " << expected_noise << std::endl;
      probabilities[k] = 0.0;
    }
  }

  return probabilities;
}

template <typename T>
Eigen::Matrix<T, 6, 1> SolveWithSnrProbabilities(
    const Eigen::Matrix<T, 6, 6>& U,
    const Eigen::Matrix<T, 6, 1>& eigenvalues,
    const Eigen::Matrix<T, 6, 1>& rhs,
    const Eigen::Matrix<T, 6, 1>& snr_probabilities) {
  typedef typename Eigen::Matrix<T, 6, 1> Vector6;

  Vector6 d_psinv = Vector6::Zero();

  for (size_t i = 0; i < 6; i++) {
    const T eigenvalue = eigenvalues[i];
    const T p = snr_probabilities[i];
    d_psinv[i] = p / eigenvalue;
  }

  Vector6 perturbation = U * d_psinv.asDiagonal() * U.transpose() * rhs;

  return perturbation;
}

template <size_t N, typename T>
auto EstimateNormal(const Eigen::Matrix<T, 3, N>& points, const T& stDevPoint, const bool& robust) {
  using Vector3 = Eigen::Matrix<T, 3, 1>;
  using Matrix3 = Eigen::Matrix<T, 3, 3>;

  Vector3 mean = Vector3::Zero();
  Matrix3 covariance = Matrix3::Zero();

  for (size_t i = 0; i < N; i++) {
    mean += points.col(i);
    covariance += points.col(i) * points.col(i).transpose();
  }
  mean /= N;
  covariance /= N;
  covariance -= mean * mean.transpose();

  Eigen::SelfAdjointEigenSolver<Matrix3> solver(covariance);
  Vector3 eigenvalues = solver.eigenvalues();
  Matrix3 eigenvectors = solver.eigenvectors();

  const Vector3 normal = eigenvectors.col(0);

  T mid_eigenvalue = eigenvalues(1);
  T max_eigenvalue = eigenvalues(2);

  if (robust) {
    mid_eigenvalue = std::max(mid_eigenvalue - stDevPoint * stDevPoint, 1e-7);
    max_eigenvalue = std::max(max_eigenvalue - stDevPoint * stDevPoint, 1e-7);
  }

  const T variance = stDevPoint * stDevPoint * (1 / T(N)) * (1 / mid_eigenvalue);
  const T distance_to_origin = normal.transpose() * mean;

  const Matrix3 covariance_of_normal = stDevPoint * stDevPoint * (1 / T(N)) * eigenvectors * Vector3(0, 1 / mid_eigenvalue, 1 / max_eigenvalue).asDiagonal() * eigenvectors.transpose();

  return std::make_tuple(normal, variance, distance_to_origin, covariance_of_normal);
}

}  // namespace degeneracy