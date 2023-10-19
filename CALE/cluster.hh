#ifndef CLUSTER_HXX
#define CLUSTER_HXX

#include <utility>
#include <vector>
#include <unordered_map>

class Cluster
{
  public:
    /**
     * @brief Constructor for a clustering algorithm.
     *
     * @param sigma
     * @param cutoff
     * @param n_rounds
     */
    Cluster(const double& sigma, const double& cutoff, const int& n_rounds);
     
    /**
     * @brief Chose a new set of particles to start clustering.   
     *
     * @param labels labels of the particles to be clustered
     * @param energies energies of the particles to be clustered
     * @param pts transverse momenta of the particles to be clustered
     * @param rapidites rapidities of the particles to be clustered
     * @param phis azimuthal angles of the particles to be clustered
     */
    void SetInputs(std::vector<int> labels,
                   std::vector<double> energies,
                   std::vector<double> pts,
                   std::vector<double> rapidites,
                   std::vector<double> phis);

    /**
     * @brief Return the labels of the current pseudojets.
     * @return vector of labels
     */
    const std::vector<int>& GetLabels() const;
    /**
     * @brief Return the energies of the current pseudojets.
     * @return vector of energies
     */
    const std::vector<double>& GetEnergies() const;
    /**
     * @brief Return the energies of the current pseudojets.
     * @param labels labels of the pseudojets to return
     * @return vector of energies
     */
    const std::vector<double> GetEnergies(const std::vector<int>& labels) const;
    /**
     * @brief Return the transverse momenta of the current pseudojets.
     * @return vector of transverse momenta
     */
    const std::vector<double>& GetPts() const;
    /**
     * @brief Return the transverse momenta of the current pseudojets.
     * @param labels labels of the pseudojets to return
     * @return vector of transverse momenta
     */
    const std::vector<double> GetPts(const std::vector<int>& labels) const;
    /**
     * @brief Return the rapidities of the current pseudojets.
     * @return vector of rapidities
     */
    const std::vector<double>& GetRapidites() const;
    /**
     * @brief Return the rapidities of the current pseudojets.
     * @param labels labels of the pseudojets to return
     * @return vector of rapidities
     */
    const std::vector<double> GetRapidites(const std::vector<int>& labels) const;
    /**
     * @brief Return the azimuthal angles of the current pseudojets.
     * @return vector of azimuthal angles
     */
    const std::vector<double>& GetPhis() const;
    /**
     * @brief Return the azimuthal angles of the current pseudojets.
     * @param labels labels of the pseudojets to return
     * @return vector of azimuthal angles
     */
    const std::vector<double> GetPhis(const std::vector<int>& labels) const;


    /**
     * @brief Say which particles form the next jet.
     * Returned as a list of labels.
     * @return vector of labels
     */
    std::vector<int> GetNextMerge() const;
    /**
     * @brief Merge a set of pesudojets into a complete jet.
     * Calculate the kinematics of the new jet internally, and
     * give it the next free label.
     *
     * @param labels labels of the pseudojets to merge
     */
    void DoMerge(std::vector<int> labels);
    /**
     * @brief Merge a set of pesudojets into a complete jet.
     * Merge the two jets, and assign the combined jet externally
     * provided label and kinematics.
     *
     * @param labels labels of the pseudojets to merge
     * @param label_new label of the new jet
     * @param energy energy of the new jet
     * @param pt transverse momentum of the new jet
     * @param rapidity rapidity of the new jet
     * @param phi azimuthal angle of the new jet
     */
    void DoMerge(std::vector<int> labels, int label_new,
                 double energy, double pt, double rapidity, double phi);

    /**
     * @brief Automatically complete the merging.
     **/
    void DoAllMerges();

    /**
     * @brief Check if all jets are complete.
     * @return true if all jets are complete
     **/
    bool IsFinished() const;

    /**
     * @brief Get a list of currently existing jets.
     *
     * Mid cluster will return the pseudojets.
     * Order matches order of GetJetConstituents, but
     * otherwise is arbitrary.
     *
     * @return each item in the vector is the label or one completed jet.
     **/
    const std::vector<int>& GetJets() const;
    /**
     * @brief Get a list of particles in each jet.
     *
     * Mid cluster will return the pseudojets.
     * Order matches order of GetJets, but otherwise is arbitrary.
     *
     * @return each item in the vector is all the labels of particles in one jet
     **/
    std::vector<std::vector<int>> GetJetConstituents() const;


    /**
     * @brief Get the next free label.
     * @return next free label
     **/
    int GetNextFreeLabel();

  private:
    /**
     * The clustering parameters
     **/
    double m_sigma;
    int m_n_rounds;
    double m_cutoff;

    /**
     * @brief A numeric label > 0 for each jet.
     **/
    std::vector<int> m_labels;
    /**
     * @brief The energy of each jet.
     **/
    std::vector<double> m_energies;
    /**
     * @brief The pt of each jet.
     **/
    std::vector<double> m_pts;
    /**
     * @brief The rapidity of each jet.
     **/
    std::vector<double> m_rapidites;
    /**
     * @brief The phi of each jet.
     **/
    std::vector<double> m_phis;

    /**
     * @brief Alternative storage of the kinematics of each jet.
     **/
    std::vector<double> m_pxs;
    std::vector<double> m_pys;
    std::vector<double> m_pzs;

    /**
     * @brief a matrix of the distances squared between particles squared.
     * By default, in the cambridge-aachen metric.
     * Used to decide seed order and also to calculate the Laplacian.
     **/
    std::vector<std::vector<double>> m_distances2;

    /**
     * @brief The calculated Laplacian of the event.
     **/
    std::vector<std::vector<double>> m_laplacian;
    /**
     * @brief The coefficients of the chebysheve polynomials.
     * This is the same for all clusterings, so it is static.
     * 50 coefficients are calculated, which is quite arbitary.
     **/
    static std::vector<double> s_chebyshev_coefficients;
    /**
     * @brief the intercal for approximation.
     * This is the same for all clusterings, so it is static.
     **/
    static std::pair<double, double> s_interval;

    /**
     * @brief the maximum number of jets the clustering could yield
     **/
    int m_max_jets;

    /**
     * @brief Indicates if a jet is avaliable for further joining.
     **/
    std::vector<bool> m_avaliable;
    /**
     * @brief Notes the index of each label in the internal list of labels.
     **/
    std::unordered_map<int, int> m_label_to_index;
    /**
     * @brief One value higher than the curant highest label value.
     **/
    int m_next_free_label = 0;
    /**
     * @brief Indicates if a jet has been merged with the beam.
     **/
    std::vector<bool> m_finished;
    /**
     * @brief The label of the parent jet for each jet.
     * Jets with no parents get a label of -1.
     **/
    std::vector<int> m_parent_labels;
    /**
     * @brief The label of each completed jet.
     **/
    std::vector<int> m_completed_labels;
    /**
     * @brief The labels of the constituents of each completed jet.
     **/
    std::vector<std::vector<int>> m_completed_constituents;

    /**
     * @brief An object to form a junk jet.
     * This is chosen as the next object that is avaliable.
     * The assumption is that all real jets have been formed.
     * If there are no available objects returns an empty vector.
     * @return vector of the label of the object to be used as a junk jet
     **/
    std::vector<int> GetSingleAvaliable() const;

    /**
     * @brief Calculate the kinematics of merging multiple jets.
     * @param labels labels of the jets to merge
     * @return vector of kinematics, energy, pt, rapidity, phi
     **/
    std::vector<double> GetMergedKinematics(std::vector<int> labels) const;

    /**
     * @brief Get the index of the object to be used as a seed.
     * If there are no seeds avaliable, returns -1.
     * @param start_seed_idx index of the seed to start from,
     *                      normally, the number of jets is a good choice
     *                      but sometimes nothing is captured by a wavlet, so the next
     *                      seed should be used.
     * @return internal index of the seed, or -1 if none avaliable
     **/
    int GetSeed(int start_seed_idx) const;

    /**
     * @brief Ordered list of indices to use as seeds
     **/
    std::vector<int> m_seed_indices;
    /**
     * @brief Total number of seeds.
     **/
    int m_n_seed_indices;
};

#endif // CLUSTER_HXX
