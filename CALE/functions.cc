#include "functions.hh"
#include <math.h>

#include <iostream>
#define MSG(x) std::cout << __FILE__ << ">" << __LINE__ << "; " << x << std::endl;
#define PRINT_VECTOR(x) std::cout << __FILE__ << ">" << __LINE__ << "; " << #x << std::endl; \
  for (auto element : x){ \
    std::cout << element << " "; \
  } \
  std::cout << std::endl;
#define PRINT_MATRIX(x) std::cout << __FILE__ << ">" << __LINE__ << "; " << #x << std::endl; \
  for (auto row : x){ \
    for (auto element : row){ \
      std::cout << element << " "; \
    } \
    std::cout << std::endl; \
  }

void Functions::RescaleLaplacian(std::vector<std::vector<double>>& laplacien){
  double scale_factor = 2./2.;
  int size = laplacien.size();
  for (int row=0; row<size; row++){
    for (int col=0; col<size; col++){
      laplacien[row][col] *= scale_factor;
    };
    laplacien[row][row] -= 1.;
  };
};

/**
 * @brief Compute Chebyshev coefficients of a kernal.
 * The kernal is a function we want to approximate in the interval
 * specified.
 * The coefficientes returned are the values we would need to multiply the
 * chebyshev series by to approximate this function.
 * @param kernal The kernal to compute the coefficients of.
 * @param max_coefficients The number of coefficients to compute.
 * @param grid_order The order of the grid to use. If -1, then max_coefficients+1 is used.
 * @param approx_interval_min lower bound of the interval of approximation.
 * @param approx_interval_max upper bound of the interval of approximation.
 * @return The Chebyshev coefficients of the kernal.
 **/
std::vector<double> Functions::ChebyshevCoefficients(double (*kernal)(const double&),
                                          const int& max_coefficients, const int& grid_order,
                                          const double& approx_interval_min, const double& approx_interval_max){
  /**
   * Formula;
   * Let f(x) be a smooth function between -1 and 1
   * Let T_n(x) be the nth chebyshev polynomial
   * Where T_0(x) = 1
   *       T_1(x) = x
   *       T_n(x) = 2xT_{n-1}(x) - T_{n-2}(x)
   * We seek the coefficients a_n such that
   * f(x) = sum_{n=0}^{m} a_n T_n(x)
   *
   * They are approximated by
   * a_n = ((2 - delta(0, n))/N) * sum_{k=0}^{N-1} f(cos(pi*(k+0.5)/N)) * T_n(cos(pi*(k+0.5)/N))
   * with a perfect approximation when N -> infinity
   **/
  double _grid_order = grid_order;
  if (grid_order == -1){
    _grid_order = max_coefficients + 1;
  };

  double scale_interval = (approx_interval_max - approx_interval_min)/2.;
  double shift_interval = (approx_interval_max + approx_interval_min)/2.;
  double grid_value, kernal_value;
  std::vector<double> grid;
  std::vector<double> kernal_values;
  double pi_over_grid_order = M_PI/_grid_order;
  for (float i=1; i<_grid_order+1; i++){
    grid_value = (i - 0.5) * pi_over_grid_order;
    grid.push_back(grid_value);
    kernal_value = kernal(std::cos(grid_value)/scale_interval+shift_interval);
    kernal_values.push_back(kernal_value);
  };
  double coefficient_value;
  std::vector<double> chebyshev_coefficients;
  double two_over_grid_order = 2./_grid_order;
  for (int i=0; i<max_coefficients+1; i++){
    coefficient_value = 0.;
    for (int j=0; j<_grid_order; j++){
      coefficient_value += kernal_values[j] * std::cos(grid[j] * i);
    };
    chebyshev_coefficients.push_back(two_over_grid_order * coefficient_value);
  };
  chebyshev_coefficients[0] /= 2.;

  return chebyshev_coefficients;
};
    
std::vector<double> Functions::ChebyshevCoefficients(const int& max_coefficients, const int& grid_order,
                                                     const double& approx_interval_min,
                                                     const double& approx_interval_max){
  return ChebyshevCoefficients([](const double& x){return std::exp(-x);}, max_coefficients, grid_order,
                               approx_interval_min, approx_interval_max);
};


std::vector<double> Functions::VectorAddition(const std::vector<double> vector1, const std::vector<double> vector2){
  int size = vector1.size();
  std::vector<double> result(size, 0.);
  for (int i=0; i<size; i++){
    result[i] = vector1[i] + vector2[i];
  };
  return result;
};

std::vector<double> Functions::VectorAddition(const double& factor1, const std::vector<double> vector1,
                                                     const double& factor2, const std::vector<double> vector2){
  int size = vector1.size();
  std::vector<double> result(size, 0.);
  for (int i=0; i<size; i++){
    result[i] = factor1*vector1[i] + factor2*vector2[i];
  };
  return result;
};


void Functions::VectorAdditionInPlace(std::vector<double>& vector1, const std::vector<double>& vector2){
  int size = vector1.size();
  for (int i=0; i<size; i++){
    vector1[i] += vector2[i];
  };
};

void Functions::VectorAdditionInPlace(const double& factor1, std::vector<double>& vector1,
                                             const double& factor2, const std::vector<double>& vector2){
  int size = vector1.size();
  for (int i=0; i<size; i++){
    vector1[i] *= factor1;
    vector1[i] += factor2 * vector2[i];
  };
};


std::vector<double> Functions::MatrixDotVector(const std::vector<std::vector<double>>& matrix,
                                                      const std::vector<double>& vector){
  int size = matrix.size();
  std::vector<double> result(size, 0.);
  for (int row=0; row<size; row++){
    for (int col=0; col<size; col++){
      result[row] += matrix[row][col] * vector[col];
    };
  };
  return result;
};

std::vector<double> Functions::LaplacianWavelet(const std::vector<std::vector<double>> &laplacian, const std::vector<double> &chebyshev_coefficients,
                                     const int& center_idx, const std::pair<double, double>& interval) {
  int n_rows = laplacian.size();
  // An all zero laplacian is a special case
  bool all_zero = true;
  for (int row = 0; row < n_rows; row++) {
    for (int col = 0; col < n_rows; col++) {
      if (laplacian[row][col] != 0.) {
        all_zero = false;
        break;
      };
    };
    if (!all_zero) break;
  };
  if (all_zero) {
    std::vector<double> result(n_rows, 0.);
    return result;
  };

  // Note, while this takes the place of make_wavelets, it's mostly cheby_op, as the starting point is the chebyshev coefficients.
  // basically equation 66 of https://inria.hal.science/hal-01943589/document
  int n_coeffients = chebyshev_coefficients.size();
  double half_interval = (interval.second - interval.first)/2.;
  double inverse_half_interval = 1./half_interval;
  double inverse_interval = 2./half_interval;
  double center = (interval.second + interval.first)/2.;
  // We will need two vectors for taking repeated dot products onto the laplacian
  // Last iteration
  std::vector<double> fourier_transform_old(n_rows, 0.);
  fourier_transform_old[center_idx] = 1.;
  // Current iteration
  std::vector<double> fourier_transform_new = VectorAddition(
      inverse_half_interval, MatrixDotVector(laplacian, fourier_transform_old),
      -center*inverse_half_interval, fourier_transform_old);
  // Placeholder for swapping them
  std::vector<double> fourier_transform_placeholder(n_rows, 0.);

  // We also need a place to store the growing sum
  std::vector<double> results = VectorAddition(0.5*chebyshev_coefficients[0], fourier_transform_old,
                                               chebyshev_coefficients[1], fourier_transform_new);

  for (int k=2; k < n_coeffients; k++){
    // update the old fourier transform to make a new one
    VectorAdditionInPlace(-1, fourier_transform_old, inverse_interval,
                          VectorAddition(1, MatrixDotVector(laplacian, fourier_transform_new),
                                         -center, fourier_transform_new));
    // switch the old and new fourier transforms
    fourier_transform_placeholder = fourier_transform_old;
    fourier_transform_old = fourier_transform_new;
    fourier_transform_new = fourier_transform_placeholder;
    // update the results
    VectorAdditionInPlace(1, results, chebyshev_coefficients[k], fourier_transform_new);
  };

  return results;
}


double Functions::AngularDistance(const double& phi1, const double& phi2){
  // taking the differnce between two angles requires careful treatment
  double delta_phi = std::abs(phi1 - phi2);
  double two_pi = 2*M_PI;
  return std::min(fmod(delta_phi, two_pi), fmod((two_pi - delta_phi), two_pi));
};


double Functions::CambridgeAachenDistance2(const double& rapidity1, const double& phi1, 
                                           const double& rapidity2, const double& phi2){
  double delta_rapidity = std::abs(rapidity1 - rapidity2);
  // taking the differnce between two angles requires careful treatment
  double delta_phi = Functions::AngularDistance(phi1, phi2);
  return delta_rapidity*delta_rapidity + delta_phi*delta_phi;
};

double Functions::GeneralisedKtDistance2(const double& pt1, const double& rapidity1, const double& phi1, 
                                         const double& pt2, const double& rapidity2, const double& phi2,
                                         const double& exponent){
  double ca_distance2 = CambridgeAachenDistance2(rapidity1, phi1, rapidity2, phi2);
  double twice_exponent = 2*exponent;
  return std::min(std::pow(pt1, twice_exponent), std::pow(pt2, twice_exponent)) * ca_distance2;
};

std::vector<std::vector<double>> Functions::GeneralisedKtDistance2Matrix(const std::vector<double>& pts,
                                                             const std::vector<double>& rapidities,
                                                             const std::vector<double>& phis,
                                                             const double& exponent){
  int n_particles = pts.size();
  std::vector<std::vector<double>> result(n_particles, std::vector<double>(n_particles, 0.));
  for (int row=0; row<n_particles; row++){
    for (int col=0; col<row; col++){
      result[row][col] = GeneralisedKtDistance2(pts[row], rapidities[row], phis[row],
                                                pts[col], rapidities[col], phis[col],
                                                exponent);
      result[col][row] = result[row][col];
    };
  };
  return result;
};

enum JetMetrics {cambridge_aachen, kt, antikt};
/**
 * @brief Distance between two particles in a named jet metric.
 * @param pt1 The transverse momentum of the first particle.
 * @param rapidity1 The rapidity of the first particle.
 * @param phi1 The azimuthal angle of the first particle.
 * @param pt2 The transverse momentum of the second particle.
 * @param rapidity2 The rapidity of the second particle.
 * @param phi2 The azimuthal angle of the second particle.
 * @param the enum value of the metric.
 * @return The distance squared between the two particles.
 **/
double Functions::NamedDistance2(const double& pt1, const double& rapidity1, const double& phi1, 
                                const double& pt2, const double& rapidity2, const double& phi2,
                                const JetMetrics& metric){
  switch (metric) {
    case cambridge_aachen:
      return CambridgeAachenDistance2(rapidity1, phi1, rapidity2, phi2);
    case kt:
      return GeneralisedKtDistance2(pt1, rapidity1, phi1, pt2, rapidity2, phi2, 1.);
    case antikt:
      return GeneralisedKtDistance2(pt1, rapidity1, phi1, pt2, rapidity2, phi2, -1.);
    default:
      throw std::invalid_argument("Invalid jet metric");
  };

};

/**
 * @brief Distance matrix between a set of particles in a named jet metric.
 * @param pts The transverse momenta of the particles.
 * @param rapidities The rapidities of the particles.
 * @param phis The azimuthal angles of the particles.
 * @param metric The enum value of the metric.
 * @return The matrix of squared distances.
 **/
std::vector<std::vector<double>> Functions::NamedDistance2Matrix(
    const std::vector<double>& pts, const std::vector<double>& rapidities, const std::vector<double>& phis,
    const JetMetrics& metric){
  switch (metric) {
    case cambridge_aachen:
      return GeneralisedKtDistance2Matrix(pts, rapidities, phis, 0.);
    case kt:
      return GeneralisedKtDistance2Matrix(pts, rapidities, phis, 1.);
    case antikt:
      return GeneralisedKtDistance2Matrix(pts, rapidities, phis, -1.);
    default:
      throw std::invalid_argument("Invalid jet metric");
  };
};


std::vector<std::vector<double>> Functions::Affinities(
        const std::vector<std::vector<double>>& distances2,
        const double& sigma){
  std::vector<std::vector<double>> affinities;
  int n_particles = distances2.size();
  // Lower triangle
  for (int row_n = 0; row_n < n_particles; row_n++){
    std::vector<double> row;
    // Calculate up to the diagonal
    for (int col_n = 0; col_n < row_n; col_n++){
      row.push_back(std::exp(-distances2[row_n][col_n]/sigma));
    };
    // The diagnal is zero
    row.push_back(0.);
    affinities.push_back(row);
  };
  // Fill in the other half of the triangle by copying
  for (int row_n = 0; row_n < n_particles; row_n++){
    for (int col_n = row_n+1; col_n < n_particles; col_n++){
      affinities[row_n].push_back(affinities[col_n][row_n]);
    };
  };
  return affinities;
};

std::vector<std::vector<double>> Functions::Laplacian(
        const std::vector<std::vector<double>>& distances2,
        const double& sigma, const bool& normalised){
  int n_particles = distances2.size();
  std::vector<std::vector<double>> laplacien = Functions::Affinities(distances2, sigma);
  // Unnormalised Laplacien
  for (int row_n = 0; row_n < n_particles; row_n++){
    double sum = 0.;
    for (int col_n = 0; col_n < n_particles; col_n++){
      sum += laplacien[row_n][col_n];
      laplacien[row_n][col_n] *= -1;
    };
    laplacien[row_n][row_n] = sum;
  };


  if (normalised){
    // Normalised Laplacien
    double diagonal, inv_sqrt_diagonal;
    for (int row_n = 0; row_n < n_particles; row_n++){
      double diagonal = laplacien[row_n][row_n];
      if (diagonal == 0.){
        inv_sqrt_diagonal = 0.;
      } else {
        inv_sqrt_diagonal = 1./std::sqrt(diagonal);
      };
      for (int col_n = 0; col_n < n_particles; col_n++){
        laplacien[row_n][col_n] *= inv_sqrt_diagonal;
        laplacien[col_n][row_n] *= inv_sqrt_diagonal;
      };
    };
  };

  return laplacien;
};

std::vector<double> Functions::PxPyPz(const double& energy, const double& pt,
                                      const double& rapidity, const double& phi){
  double px = pt*std::cos(phi);
  double py = pt*std::sin(phi);
  double rapidity_factor = std::exp(2*rapidity);
  double pz = energy*(rapidity_factor-1)/(1+rapidity_factor);
  return {px, py, pz};
};

std::vector<double> Functions::PtRapPhi(const double& energy, const double& px,
                                        const double& py, const double& pz){
  double pt = std::sqrt(px*px + py*py);
  double phi = std::atan2(py, px);
  double rapidity = std::copysign(0.5*std::log((energy+pz)/(energy-pz)), pz);
  return {pt, rapidity, phi};
};




