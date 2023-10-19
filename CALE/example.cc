// To run this example, use the following command:
//
//   ./example < THE_INPUT_FILE [replace this with the file you want used!]
//
//----------------------------------------------------------------------
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

#include <iostream>
#include <sstream>

#include "fastjet/PseudoJet.hh"
#include <sstream>
#include "CALE.hh" // In external code, this should be fastjet/contrib/CALE.hh

using namespace std;
using namespace fastjet;
using namespace contrib;

void print_jets (const fastjet::ClusterSequence &,
                 const vector<fastjet::PseudoJet> &);


// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

//----------------------------------------------------------------------
int main(){

  //----------------------------------------------------------
  // read in input particles
  vector<PseudoJet> event;
  read_event(event);
  cout << "# read an event with " << event.size() << " particles" << endl;

  //----------------------------------------------------------
  // illustrate how this CALE contrib works
  // defining parameters
  double sigma = 0.1;
  double cutoff = 0.0;
  int n_rounds = 15;
  
  CALEPlugin jet_pluginCALE(sigma, cutoff, n_rounds);
  fastjet::JetDefinition jet_def(&jet_pluginCALE);
  fastjet::ClusterSequence clust_seq(event, jet_def);
  
  // tell the user what was done
  cout << "# Ran " << jet_def.description() << endl;

  // extract the inclusive jets with pt > 5 GeV
  double ptmin = 5.0;
  vector<fastjet::PseudoJet> inclusive_jets = clust_seq.inclusive_jets(ptmin);
  
  // print them out
  cout << "Printing inclusive jets with pt > "<< ptmin <<" GeV\n";
  cout << "---------------------------------------\n";
  print_jets(clust_seq, inclusive_jets);
  cout << endl;
  
  return 0;
}

// read in input particles
void read_event(vector<PseudoJet> &event){  
  string line;
  while (getline(cin, line)) {
    istringstream linestream(line);
    // take substrings to avoid problems when there is extra "pollution"
    // characters (e.g. line-feed).
    if (line.substr(0,4) == "#END") {return;}
    if (line.substr(0,1) == "#") {continue;}
    double px,py,pz,E;
    linestream >> px >> py >> pz >> E;
    PseudoJet particle(px,py,pz,E);

    // push event onto back of full_event vector
    event.push_back(particle);
  }
}

//----------------------------------------------------------------------
/// a function that pretty prints a list of jets
void print_jets (const fastjet::ClusterSequence & clust_seq,
                 const vector<fastjet::PseudoJet> & jets) {
  
  // sort jets into increasing pt
  vector<fastjet::PseudoJet> sorted_jets = sorted_by_pt(jets);
  
  // label the columns
  printf("%5s %10s %10s %10s %10s %10s %10s\n","jet #", "rapidity",
         "phi", "pt","m","e", "n constituents");
  
  // print out the details for each jet
  for (unsigned int i = 0; i < sorted_jets.size(); i++) {
    int n_constituents = clust_seq.constituents(sorted_jets[i]).size();
    printf("%5u %10.3f %10.3f %10.3f %10.3f %10.3f %8u\n",
           i, sorted_jets[i].rap(), sorted_jets[i].phi(),
           sorted_jets[i].perp(),sorted_jets[i].m(),sorted_jets[i].e(), n_constituents);
  }
}
