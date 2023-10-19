#include "cluster.hh"
#include "functions.hh"
#include <cassert>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <stdexcept>

#include <iostream>
#define MSG(x) std::cout << __FILE__ << ">" << __LINE__ << "; " << x << std::endl;
#define PRINT_VECTOR(x) std::cout << __FILE__ << ">" << __LINE__ << "; " << #x << std::endl; \
  for (auto element : x){ \
    std::cout << element << " "; \
  } \
  std::cout << std::endl << std::endl;
#define PRINT_MATRIX(x) std::cout << __FILE__ << ">" << __LINE__ << "; " << #x << std::endl; \
  for (auto row : x){ \
    for (auto element : row){ \
      std::cout << element << " "; \
    } \
    std::cout << std::endl; \
  } \
  std::cout << std::endl;

Cluster::Cluster(const double& sigma, const double& cutoff, const int& n_rounds) :
  m_sigma(sigma), m_cutoff(cutoff), m_n_rounds(n_rounds) {};


// Calculate the static chebyshev coefficients.
std::vector<double> Cluster::s_chebyshev_coefficients = Functions::ChebyshevCoefficients(50);
// The approximation interval is also static
// Becuase of normalisation, the max eigenvalue is 2., so the interval is [0, 2]
// The chebyshev polynomials are defined over -1 to 1, so that range needs to be shifted.
// See https://inria.hal.science/hal-01943589/document page 23
std::pair<double, double> Cluster::s_interval = std::make_pair(0., 2.);


void Cluster::SetInputs(std::vector<int> labels,
                        std::vector<double> energies,
                        std::vector<double> pts,
                        std::vector<double> rapidites,
                        std::vector<double> phis){
  m_labels = labels;
  m_energies = energies;
  m_pts = pts;
  m_rapidites = rapidites;
  m_phis = phis;
  // Check all given properties are the same length
  int n_labels = labels.size();
  assert(energies.size() == n_labels);
  assert(pts.size() == n_labels);
  assert(rapidites.size() == n_labels);
  assert(phis.size() == n_labels);
  for (int i=0; i<n_labels; i++){
    // Create a map that makes for quick lookups of the index of a label
    m_label_to_index[labels[i]] = i;
    // Update next free label
    if (labels[i] >= m_next_free_label){
      m_next_free_label = labels[i] + 1;
    }
    // Calculate the px, py, pz of each particle
    std::vector<double> pxpypx = Functions::PxPyPz(energies[i], pts[i], rapidites[i], phis[i]);
    m_pxs.push_back(pxpypx[0]);
    m_pys.push_back(pxpypx[1]);
    m_pzs.push_back(pxpypx[2]);
  }
  // All particles become pseudojets and are avaliable for merging
  m_avaliable = std::vector<bool>(n_labels, true);
  // No jets are finished
  m_finished = std::vector<bool>(n_labels, false);
  // No pseudojets have parents yet
  m_parent_labels = std::vector<int>(n_labels, -1);
  // Calculate the laplacian
  m_distances2 = Functions::NamedDistance2Matrix(m_pts, m_rapidites, m_phis,
                                                 Functions::cambridge_aachen);
  m_laplacian = Functions::Laplacian(m_distances2, m_sigma, true);
  // Check the max number of jets
  m_max_jets = m_n_rounds < n_labels ? m_n_rounds : n_labels;
  // Decide on the seed order
  // This is done with akt, not cambridge-aachen.
  std::vector<std::vector<double>> antikt_distances2 = Functions::NamedDistance2Matrix(
      m_pts, m_rapidites, m_phis, Functions::antikt);
  std::vector<double> summed_distances = std::vector<double>(n_labels, 0.);
  for (int i=0; i<n_labels; i++){
    for (int j=0; j<n_labels; j++){
      summed_distances[i] += std::sqrt(antikt_distances2[i][j]);
    }
  }
  m_seed_indices = std::vector<int>(n_labels);
  m_n_seed_indices = m_seed_indices.size();

  std::iota(m_seed_indices.begin(), m_seed_indices.end(), 0);
  std::sort(m_seed_indices.begin(), m_seed_indices.end(),
            [&summed_distances](int i1, int i2){return summed_distances[i1] < summed_distances[i2];});
};

const std::vector<int>& Cluster::GetLabels() const {
  return m_labels;
};

const std::vector<double>& Cluster::GetEnergies() const {
  return m_energies;
};

const std::vector<double> Cluster::GetEnergies(const std::vector<int>& labels) const {
  std::vector<double> energies(labels.size());
  for (int i=0; i<labels.size(); i++){
    energies[i] = m_energies[m_label_to_index.at(labels[i])];
  }
  return energies;
};

const std::vector<double>& Cluster::GetPts() const {
  return m_pts;
};

const std::vector<double> Cluster::GetPts(const std::vector<int>& labels) const {
  std::vector<double> pts(labels.size());
  for (int i=0; i<labels.size(); i++){
    pts[i] = m_pts[m_label_to_index.at(labels[i])];
  }
  return pts;
};

const std::vector<double>& Cluster::GetRapidites() const {
  return m_rapidites;
};

const std::vector<double> Cluster::GetRapidites(const std::vector<int>& labels) const {
  std::vector<double> rapidites(labels.size());
  for (int i=0; i<labels.size(); i++){
    rapidites[i] = m_rapidites[m_label_to_index.at(labels[i])];
  }
  return rapidites;
};

const std::vector<double>& Cluster::GetPhis() const {
  return m_phis;
};

const std::vector<double> Cluster::GetPhis(const std::vector<int>& labels) const {
  std::vector<double> phis(labels.size());
  for (int i=0; i<labels.size(); i++){
    phis[i] = m_phis[m_label_to_index.at(labels[i])];
  }
  return phis;
};

int Cluster::GetSeed(int start_seed_idx) const {
  // The seed idx must be at least the number of completed jets in
  int idx = start_seed_idx;
  if (idx >= m_n_seed_indices){
    return -1;
  }
  while (m_avaliable[m_seed_indices[idx]] == false){
    idx++;
    if (idx >= m_n_seed_indices){
      return -1;
    }
  }  
  return m_seed_indices[idx];
};

std::vector<int> Cluster::GetSingleAvaliable() const{
  for (int i=0; i<m_avaliable.size(); i++){
    if (m_avaliable[i] == true){
      return {m_labels[i]};
    }
  }
  return {};
};


std::vector<int> Cluster::GetNextMerge() const {
  std::vector<int> labels;
  int start_seed_idx = m_completed_labels.size();
  int seed;  // Get inside the loop so we update if a seed doesn't work
  while (labels.size() == 0){
    seed = this->GetSeed(start_seed_idx);
    if (seed == -1){
      // We have done all the real clustering,
      // just return the first avaliable as a junk jet.
      return GetSingleAvaliable();
    }
    // Decide what is close to the seed
    std::vector<double> wavelets = Functions::LaplacianWavelet(m_laplacian, 
                                                               Cluster::s_chebyshev_coefficients,
                                                               seed, 
                                                               Cluster::s_interval);
    double max_wavelet = *std::max_element(wavelets.begin(), wavelets.end());
    double min_wavelet = *std::min_element(wavelets.begin(), wavelets.end());
    double shifted_threshold = (max_wavelet - min_wavelet)*(m_cutoff + 1.)/2. + min_wavelet;
    for (int i=0; i<wavelets.size(); i++){
      // Don't bother with things that aren't available
      if (!m_avaliable[i]){
        continue;
      }
      // Scale the wavelets from -1 to 1
      if (wavelets[i] < shifted_threshold){
        labels.push_back(m_labels[i]);
      }
    }
    // If no labels were found, we need to increase the seed index
    start_seed_idx++;
  }
  return labels;
};


void Cluster::DoMerge(std::vector<int> labels){
  int label_new = GetNextFreeLabel();
  std::vector<double> kinematics = GetMergedKinematics(labels);
  DoMerge(labels, label_new, kinematics[0], kinematics[1], kinematics[2], kinematics[3]);
};


void Cluster::DoMerge(std::vector<int> labels, int label_new,
                      double energy, double pt, double rapidity, double phi){
  // First log the jet completion
  m_completed_labels.push_back(label_new);
  m_next_free_label = label_new + 1;
  m_completed_constituents.push_back(labels);
  // Then update the listed labels
  m_labels.push_back(label_new);
  m_energies.push_back(energy);
  m_pts.push_back(pt);
  m_rapidites.push_back(rapidity);
  m_phis.push_back(phi);
  // Update the label map
  m_label_to_index[label_new] = m_labels.size()-1;
  // Get the px, py, pz of the new particle
  std::vector<double> pxpypz = Functions::PxPyPz(energy, pt, rapidity, phi);
  m_pxs.push_back(pxpypz[0]);
  m_pys.push_back(pxpypz[1]);
  m_pzs.push_back(pxpypz[2]);
  // update the availability and parantage
  for (int label : labels){
    m_avaliable[m_label_to_index[label]] = false;
    m_parent_labels[m_label_to_index[label]] = label_new;
  }
  // The final jet is not only unavailable, but also finished
  m_avaliable.push_back(false);
  m_finished.push_back(true);
  m_parent_labels.push_back(-1);
};

void Cluster::DoAllMerges(){
  while (!IsFinished()){
    std::vector <int> merge = GetNextMerge();
    DoMerge(merge);
  }
};

bool Cluster::IsFinished() const {
  if (m_completed_labels.size() == m_max_jets){
    return true;
  }
  // We do not expect all pseudojets to be finished, the ones
  // that are merged into other jets are "unfinished"
  // but we do expect that none of them are available
  return std::none_of(m_avaliable.begin(), m_avaliable.end(), [](bool v) { return v; });
};

const std::vector<int>& Cluster::GetJets() const {
  return m_completed_labels;
};


std::vector<std::vector<int>> Cluster::GetJetConstituents() const {
  return m_completed_constituents;
};


int Cluster::GetNextFreeLabel(){
  return m_next_free_label;
};


std::vector<double> Cluster::GetMergedKinematics(std::vector<int> labels) const {
  double e = 0;
  double px = 0;
  double py = 0;
  double pz = 0;
  int idx;
  for (int label : labels){
    idx = m_label_to_index.at(label); 
    e += m_energies[idx];
    px += m_pxs[idx];
    py += m_pys[idx];
    pz += m_pzs[idx];
  }
  std::vector<double> kinematics = Functions::PtRapPhi(e, px, py, pz);
  kinematics.insert(kinematics.begin(), e);
  return kinematics;
};
