// main91.cc is a part of the PYTHIA event generator.
// Copyright (C) 2020 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Keywords: analysis; root;

// This is a simple test program.
// It studies the charged multiplicity distribution at the LHC.
// Modified by Rene Brun, Axel Naumann and Bernhard Meirose
// to use ROOT for histogramming.

// Stdlib header file for input and output.
#include <iostream>
#include <algorithm>
#include <vector>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include <cmath>
#include "helpers.hh"
#include "Pythia8/FJcore.h"

// ROOT, for histogramming.
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"

using namespace Pythia8;
using namespace std;
using namespace fjcore;

bool myfunction (double a, double b) {return (a<b); }

int main(int argc, char* argv[]) {

  // Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);

  // Create Pythia instance and set it up to generate hard QCD processes
  // above pTHat = 20 GeV for pp collisions at 14 TeV.
  Pythia pythia;
  pythia.readString("Beams:eCM = 91.2");
  pythia.readString("Beams:idA = 11"); //electron
  pythia.readString("Beams:idB = -11"); //positron
  
  // produce Z bosons and interference with photon
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");

  //consider only Z (PDGID=23) decays to light quarks
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");

  // by changing the seed we get different events
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 20");
  pythia.init();

  // Create file on which histogram(s) can be saved.
  TFile* outFile = new TFile("jetHistos.root", "RECREATE");

  // Book histogram.
  TH1F *mjet = new TH1F("mjet", "Mass of leading jet in E",90, 0., 90.);
  TH1F *ejet = new TH1F("ejet", "Energy of leading jet in E", 90, 0., 90.);
  TH1F *meratio = new TH1F("meratio", "M/E ratio of leading jet in E", 30, 0, 3);
  TH1F *npar = new TH1F("npar", "Number of constituent particles in leading jet", 40, 0, 40);

  int nev = 10000;
  // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nev; iEvent++){
    if (!pythia.next()) continue;

    vector<PseudoJet> particles;

    for(int i = 0; i < pythia.event.size(); i++){

      // Consider only final state particles
      if(!pythia.event[i].isFinal()) continue;

     particles.push_back(PseudoJet(pythia.event[i]));

    }// End of particles loop

  //Cluster particles into jets
  //First generate a whole "clustering sequence" with the e+e- kt algorithm
    JetDefinition jet_def(ee_kt_algorithm);
    ClusterSequence cs(particles, jet_def);
    //then select the jets that come out of clustering with the following
    // value of ycut

    double ycut = 0.02;
    vector<PseudoJet> jets = cs.exclusive_jets_ycut(ycut);

    //d::sort(jets.begin(), jets.end());
 //std::cout << "Event number: " << iEvent <<std::endl;
//  for(unsigned i = 0; i < jets.size(); i++){
//    std::cout << "Energy: " << jets[i].e() << " pT: " << jets[i].pt() << std::endl;
//  }
  int njet = jets.size();

  mjet->Fill(jets[njet-1].m());
  ejet->Fill(jets[njet-1].e());
  meratio->Fill(jets[njet-1].m()/jets[njet-1].e());
  npar->Fill(jets[njet-1].constituents().size());

    
  }//End of event loop

  // Statistics on event generation.
  pythia.stat();

  // Show histogram. Possibility to close it.
 //ult->Draw();
  mjet->Draw();
  ejet->Draw();
  meratio->Draw();
  npar->Draw();
//std::cout << "\nDouble click on the histogram window to quit.\n";
//gPad->WaitPrimitive();

  // Save histogram on file and close file.
//mult->Write();
  mjet->Write();
  ejet->Write();
  meratio->Write();
  npar->Write();

  delete outFile;

  // Done.
  return 0;
}
