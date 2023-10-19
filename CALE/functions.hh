#ifndef FUNCTIONS_HXX
#define FUNCTIONS_HXX

#include <vector>

// Making a class of this, in case future modifications require some
// clever use of existing calculations.
class Functions
{
  public:


    /**
     * @brief Rescale eigenvalues of the laplacien to -1 to 1.
     * Operates in place.
     * @param laplacien The laplacien.
     */
    static void RescaleLaplacian(std::vector<std::vector<double>>& laplacien);

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
    static std::vector<double> ChebyshevCoefficients(double (*kernal)(const double&),
                                              const int& max_coefficients, const int& grid_order=-1,
                                              const double& approx_interval_min=-1.0, const double& approx_interval_max=1.0);
        
    /**
     * @brief Compute Chebyshev coefficients of a negative exponent kernal.
     * The coefficientes returned are the values we would need to multiply the
     * chebyshev series by to approximate this exp(-x)
     * @param max_coefficients The number of coefficients to compute.
     * @param grid_order The order of the grid to use. If -1, then max_coefficients+1 is used.
     * @param approx_interval_min lower bound of the interval of approximation.
     * @param approx_interval_max upper bound of the interval of approximation.
     * @return The Chebyshev coefficients of the kernal.
     **/
    static std::vector<double> ChebyshevCoefficients(const int& max_coefficients, const int& grid_order=-1,
                                              const double& approx_interval_min=0.0, const double& approx_interval_max=2.0);

    /**
     * @brief Sum two vectors.
     * @param vector1 The first vector.
     * @param vector2 The second vector.
     * @return vector1 + vector2.
     **/
    static std::vector<double> VectorAddition(const std::vector<double> vector1, const std::vector<double> vector2);

    /**
     * @brief Sum two vectors with rescaling.
     * @param factor1 The factor to multiply vector1 by.
     * @param vector1 The first vector.
     * @param factor2 The factor to multiply vector2 by.
     * @param vector2 The second vector.
     * @return vector1 + vector2.
     **/
    static std::vector<double> VectorAddition(const double& factor1, const std::vector<double> vector1,
                                              const double& factor2, const std::vector<double> vector2);

    /**
     * @brief Sum two vectors in place.
     * @param vector1 The first vector, which will be altered.
     * @param vector2 The second vector, which will be subtracted from the first.
     **/
    static void VectorAdditionInPlace(std::vector<double>& vector1, const std::vector<double>& vector2);

    /**
     * @brief Sum two vectors in place.
     * @param factor1 The factor to multiply vector1 by.
     * @param vector1 The first vector, which will be altered.
     * @param factor2 The factor to multiply vector2 by.
     * @param vector2 The second vector, which will be subtracted from the first.
     **/
    static void VectorAdditionInPlace(const double& factor1, std::vector<double>& vector1,
                                      const double& factor2, const std::vector<double>& vector2);

    /**
     * @brief Dot product between a matrix and a vector
     * @param matrix The matrix.
     * @param vector The vector.
     * @return The dot product.
     **/
    static std::vector<double> MatrixDotVector(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector);

    /**
     * @brief Using Chebyshev polynomials, find points with similar eigenvector components to the center_idx.
     * @param laplacian The laplacian.
     * @param chebyshev_coefficients The Chebyshev coefficients.
     * @param center_idx The index to center the approximation around.
     * @param interval The interval to approximate in.
     * @return The approximation of the eigenvalues.
     **/
    static std::vector<double> LaplacianWavelet(const std::vector<std::vector<double>> &laplacian,
                                         const std::vector<double> &chebyshev_coefficients,
                                         const int& center_idx, const std::pair<double, double>& interval);
    /**
     * @brief Distance between to angles.
     * @param phi1 The first angle.
     * @param phi2 The second angle.
     * @return The distance between the two angles.
     **/
    static double AngularDistance(const double& phi1, const double& phi2);

    /**
     * @brief Distance squared between two particles in the Cambridge-Aachen metric.
     * @param rapidity1 The rapidity of the first particle.
     * @param phi1 The azimuthal angle of the first particle.
     * @param rapidity2 The rapidity of the second particle.
     * @param phi2 The azimuthal angle of the second particle.
     * @return The distance squared between the two particles.
     **/
    static double CambridgeAachenDistance2(const double& rapidity1, const double& phi1, 
                                    const double& rapidity2, const double& phi2);

    /**
     * @brief Distance squared between two particles in the Generalised-Kt metric.
     * @param pt1 The transverse momentum of the first particle.
     * @param rapidity1 The rapidity of the first particle.
     * @param phi1 The azimuthal angle of the first particle.
     * @param pt2 The transverse momentum of the second particle.
     * @param rapidity2 The rapidity of the second particle.
     * @param phi2 The azimuthal angle of the second particle.
     * @param the exponent of the metric.
     * @return The distance squared between the two particles.
     **/
    static double GeneralisedKtDistance2(const double& pt1, const double& rapidity1, const double& phi1, 
                                 const double& pt2, const double& rapidity2, const double& phi2,
                                 const double& exponent);

    /**
     * @brief Distance squared matrix between a set of particles in the Generalised-Kt metric.
     * Note, the distance is squared, not the matrix.
     * @param pts The transverse momenta of the particles.
     * @param rapidities The rapidities of the particles.
     * @param phis The azimuthal angles of the particles.
     * @param the exponent of the metric.
     * @return The distance squared matrix between the particles.
     **/
    static std::vector<std::vector<double>> GeneralisedKtDistance2Matrix(const std::vector<double>& pts,
                                                                 const std::vector<double>& rapidities,
                                                                 const std::vector<double>& phis,
                                                                 const double& exponent);

    enum JetMetrics {cambridge_aachen, kt, antikt};
    /**
     * @brief Distance squared between two particles in a named jet metric.
     * @param pt1 The transverse momentum of the first particle.
     * @param rapidity1 The rapidity of the first particle.
     * @param phi1 The azimuthal angle of the first particle.
     * @param pt2 The transverse momentum of the second particle.
     * @param rapidity2 The rapidity of the second particle.
     * @param phi2 The azimuthal angle of the second particle.
     * @param the enum value of the metric.
     * @return The distance squared between the two particles.
     **/
    static double NamedDistance2(const double& pt1, const double& rapidity1, const double& phi1, 
                         const double& pt2, const double& rapidity2, const double& phi2,
                         const JetMetrics& metric);

    /**
     * @brief Distance squared matrix between a set of particles in a named jet metric.
     * Note, the distance is squared, not the matrix.
     * @param pts The transverse momenta of the particles.
     * @param rapidities The rapidities of the particles.
     * @param phis The azimuthal angles of the particles.
     * @param metric The enum value of the metric.
     * @return The matrix of squared distances.
     **/
    static std::vector<std::vector<double>> NamedDistance2Matrix(
        const std::vector<double>& pts, const std::vector<double>& rapidities, const std::vector<double>& phis,
        const JetMetrics& metric);

    /**
     * @brief Affinity matrix for a set of particles.
     * Assumes distance2 is symmetric.
     * @param distances2 The distances squared in the chosen metric between particles.
     * @param sigma The sigma in the denominator of the exponent of the affinity.
     * @return The matrix of affinites.
     **/
    static std::vector<std::vector<double>> Affinities(
        const std::vector<std::vector<double>>& distances2,
        const double& sigma);

    /**
     * @brief Graph laplacien for a set of particles.
     * @param distances2 The distances squared in the chosen metric between particles.
     * @param sigma The sigma in the denominator of the exponent of the affinity.
     * @param normalised Whether to normalise the laplacien.
     * @return The laplacien.
     **/
    static std::vector<std::vector<double>> Laplacian(
        const std::vector<std::vector<double>>& distances2,
        const double& sigma, const bool& normalised=true);

    /**
     * @brief From detector coordinates, calculate cartesien coordinates
     * @param energy The energy.
     * @param pt The transverse momentum.
     * @param rapidity The rapidity.
     * @param phi The azimuthal angle.
     * @return The px, py, pz.
     **/
    static std::vector<double> PxPyPz(const double& energy, const double& pt,
                               const double& rapidity, const double& phi);

    /**
     * @brief From cartesien coordinates, calculate detector coordinates
     * @param energy The energy.
     * @param px The x momentum.
     * @param py The y momentum.
     * @param pz The z momentum.
     * @return The pt, rapidity, phi.
     **/
    static std::vector<double> PtRapPhi(const double& energy, const double& px, const double& py, const double& pz);

};

#endif // FUNCTIONS_HXX
