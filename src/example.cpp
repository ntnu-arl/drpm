#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include "data.h"
#include "degeneracy.h"

Eigen::Matrix<double, 6, 6> ComputeHessian(const degeneracy::VectorVector3<double>& points, const degeneracy::VectorVector3<double>& normals, const std::vector<double>& weights) {
  const size_t nPoints = points.size();
  Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero(6, 6);
  for (size_t i = 0; i < nPoints; i++) {
    const Eigen::Vector3d point = points[i];
    const Eigen::Vector3d normal = normals[i];
    const Eigen::Vector3d pxn = point.cross(normal);
    const double w = std::sqrt(weights[i]);
    Eigen::Matrix<double, 6, 1> v;
    v.head(3) = w * pxn;
    v.tail(3) = w * normal;
    H += v * v.transpose();
  }
  return H;
}

degeneracy::VectorMatrix3<double> GetIsotropicCovariances(const size_t& N, const double stdev) {
  degeneracy::VectorMatrix3<double> covariances;
  covariances.reserve(N);
  for (size_t i = 0; i < N; i++) {
    covariances.push_back(Eigen::Matrix3d::Identity() * std::pow(stdev, 2));
  }
  return covariances;
}

int main() {
  // Points, normals and covariances must be expressed in the same frame of reference
  // For the conditioning of the Hessian, it is preferable to use the LiDAR frame (and not the world frame)
  const auto points = data::points;
  const auto normals = data::normals;
  const auto weights_squared = data::weights_squared;
  const auto normal_covariances = GetIsotropicCovariances(data::normals.size(), data::stdev_normals);

  const auto H = ComputeHessian(points, normals, weights_squared);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> eigensolver(H);

  const auto eigenvectors = eigensolver.eigenvectors();
  const auto eigenvalues = eigensolver.eigenvalues();

  Eigen::Matrix<double, 6, 6> noise_mean;
  Eigen::Matrix<double, 6, 1> noise_variance;
  const double snr_factor = 10.0;

  std::tie(noise_mean, noise_variance) = degeneracy::ComputeNoiseEstimate<double, double>(points, normals, weights_squared, normal_covariances, eigenvectors, data::stdev_points);
  Eigen::Matrix<double, 6, 1> non_degeneracy_probabilities = degeneracy::ComputeSignalToNoiseProbabilities<double>(H, noise_mean, noise_variance, eigenvectors, snr_factor);

  std::cout << "The non-degeneracy probabilities are: " << std::endl;
  std::cout << non_degeneracy_probabilities.transpose() << std::endl;

  std::cout << "For the eigenvectors of the Hessian: " << std::endl;
  std::cout << eigenvectors << std::endl;

  // // The following exemplifies how to solve the system of equations using the probabilities
  // // Dummy right hand side rhs = Jtb
  // const Eigen::Matrix<double, 6, 1> rhs = Eigen::Matrix<double, 6, 1>::Zero(6, 1);
  // const auto estimate = degeneracy::SolveWithSnrProbabilities(eigenvectors, eigenvalues, rhs, non_degeneracy_probabilities);

  return 0;
}
