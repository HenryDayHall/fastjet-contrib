// $Id$
//
// Copyright (c) -, 
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include "CALE.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{


// constructor 
CALEPlugin::CALEPlugin(const double& sigma, const double& cutoff, const int& n_rounds):
  m_sigma(sigma), m_cutoff(cutoff), m_n_rounds(n_rounds),
  m_algorithm(sigma, cutoff, n_rounds){};


std::string CALEPlugin::description() const{
  // a description includes the algorithm parameters
  std::ostringstream ostr;
  ostr << "jet algorithm is Chebyshev Approximated Laplacian Eigenvectors: sigma = " << m_sigma << ", cutoff = " << m_cutoff << ", n_rounds = " << m_n_rounds;
  return ostr.str();
};


void CALEPlugin::run_clustering(ClusterSequence &cs) const{
  m_cs = &cs;  // update the cluster sequence
  std::vector<int> labels;
  std::vector<double> energies, pts, rapidities, phis;
  for ( const auto jet : cs.jets() ){
    labels.push_back(jet.cluster_hist_index());
    energies.push_back(jet.e());
    pts.push_back(jet.pt());
    rapidities.push_back(jet.rapidity());
    phis.push_back(jet.phi());
  }
  m_algorithm.SetInputs(labels, energies, pts, rapidities, phis);
  m_algorithm.DoAllMerges();
  std::vector<std::vector<int>> jet_constituents = m_algorithm.GetJetConstituents();
  for (std::vector<int> this_jet : jet_constituents){
    // These objects are all in one jet
    MakeJet(this_jet);
  }
};

int CALEPlugin::MakeJet(const std::vector<int> &cluster_hist_indices){
  double distance = -1.;  // We don't have real distances to record.
  int n_leaves = cluster_hist_indices.size();
  int last_index_created = cluster_hist_indices[0];
  int next_index;
  for (int leaf_n=1; leaf_n<n_leaves; leaf_n++){
    m_cs.plugin_record_ij_recombination(last_index_created, cluster_hist_indices[leaf_n], distance, next_index);
    last_index_created = next_index;
  }
  m_cs.plugin_record_iB_recombination(last_index_created, distance);
  return 0;
};

} // namespace contrib

FASTJET_END_NAMESPACE
