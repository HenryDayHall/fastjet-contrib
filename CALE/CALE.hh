// $Id$
//
// Copyright (c) 2023, Sirnandan Dasmahapatra, Henry Day-Hall,
// Kieran Maguire and Stefano Moretti.
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

#ifndef __FASTJET_CONTRIB_CALE_HH__
#define __FASTJET_CONTRIB_CALE_HH__

#include <fastjet/internal/base.hh>

// My includes
#include "cpp_CALE/CALE/cluster.hxx"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//------------------------------------------------------------------------
/// \class CALEPlugin
/// <insert short description>
///
/// <lnsert long description>
class CALEPlugin : public JetDefinition::Plugin {
public:
  /// Constructor
  /// \param sigma algorithm parameter
  /// \param cutoff algorithm parameter
  /// \param n_rounds algorithm parameter
  CALEPlugin(const double& sigma, const double& cutoff, const int& n_rounds);
  // Methods required by base class.
  /// Description of CALE, inc parameters
  virtual std::string description() const;
  /// Cluster all particles, inc. noise as single particle jets.
  virtual void run_clustering(ClusterSequence &cs) const;
  /// This is a jet parameter. The fastjet library assumes you just have the one.
  /// We will consider it to be the Cutoff value, and Sigma and NRounds 
  virtual double R() const {return m_cutoff;};
  /// Exclusive jets are currently the only option going.
  /// Though, in principle, inclusive clustering would be possible
  /// with the same aproach.
  virtual bool exclusive_sequence_meaningful() const {return true};

private:

  /// Algorithm is implemented in a seperate class
  /// which we wrap here
  Cluster m_algorithm;
  /// Algorithm parameters
  double m_sigma;
  double m_cutoff;
  int m_n_rounds;
  /// The current cluster sequence
  ClusterSequence m_current;
  /// Function to merge a set of particle into a jet (given the cluster history index)
  /// \param cluster_hist_indices the indices of the leaf particles of the jet
  /// \return the cluster history index of the jet
  int MakeJet(const std::vector<int> &cluster_hist_indices);
};


} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_CALE_HH__
