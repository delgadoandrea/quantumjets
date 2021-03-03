#include "Pythia8/Pythia.h"
#include <string>
#include <sstream>

using namespace Pythia8;

int main() {

  // Generator.
  Pythia pythia;


  // Allow no substructure in e+- beams: normal for corrected LEP data.
  pythia.readString("PDF:lepton = off");
  // Process selection.
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  // Switch off all Z0 decays and then switch back on those to quarks.
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");

  // LEP1 initialization at Z0 mass.
  pythia.readString("Beams:idA =  11");
  pythia.readString("Beams:idB = -11");
  double mZ = pythia.particleData.m0(23);
  pythia.settings.parm("Beams:eCM", mZ);
  pythia.init();

  // Set up Sphericity, "Linearity", Thrust and cluster jet analyses.
  Sphericity sph;
  Sphericity lin(1.);
  Thrust thr;
  ClusterJet lund("Lund");
  ClusterJet jade("Jade");
  ClusterJet durham("Durham");

  std::string path = "/home/andrea/pythiaEvents/";

  std::string fres = path + "Pythia8Results.dat";
  ofstream filer(fres);


  // Begin event loop. Generate event. Skip if error. List first few.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;

    std::string fout = path + "Event_" + std::to_string(iEvent) + "_v2.dat";
    ofstream file(fout);


    // Find and histogram charged multiplicity.
    int nCh = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) ++nCh;

    // Calculate sphericity.
    if (sph.analyze( pythia.event )) {
      if (iEvent < 3) sph.list();
      double e1 = sph.eigenValue(1);
      double e2 = sph.eigenValue(2);
      double e3 = sph.eigenValue(3);
      if (e2 > e1 || e3 > e2) cout << "eigenvalues out of order: "
      << e1 << "  " << e2 << "  " << e3 << endl;
    }

    // Find linearized sphericity.
    if (lin.analyze( pythia.event )) {
      if (iEvent < 3) lin.list();
      double e1 = lin.eigenValue(1);
      double e2 = lin.eigenValue(2);
      double e3 = lin.eigenValue(3);
      if (e2 > e1 || e3 > e2) cout << "eigenvalues out of order: "
      << e1 << "  " << e2 << "  " << e3 << endl;
    }

    // Find and histogram thrust.
    if (thr.analyze( pythia.event )) {
      if (iEvent < 3) thr.list();
      if ( abs(thr.eventAxis(1).pAbs() - 1.) > 1e-8
        || abs(thr.eventAxis(2).pAbs() - 1.) > 1e-8
        || abs(thr.eventAxis(3).pAbs() - 1.) > 1e-8
        || abs(thr.eventAxis(1) * thr.eventAxis(2)) > 1e-8
        || abs(thr.eventAxis(1) * thr.eventAxis(3)) > 1e-8
        || abs(thr.eventAxis(2) * thr.eventAxis(3)) > 1e-8 ) {
        cout << " suspicious Thrust eigenvectors " << endl;
        thr.list();
      }
    }

    // Calculate alternate version of thrust by using linear sphericity as seed reference axis

      // Initial values and counters zero.
    double eVal1 = 0.;
    double tMax = 0.0;
    double firstVal = 0.0;
    Vec4 eVec1 = 0.;
    int nStudy = 0;
    int nFew = 0;
    vector<Vec4> pOrder;
    Vec4 pSum, aRef, pPart, pFull, pMax;
    int nIter = 0;

    // Loop over desired particles in the event.
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].isFinal()) {
        //if (select >  2 &&  pythia.event[i].isNeutral() ) continue;
        //if (select == 2 && !pythia.event[i].isVisible() ) continue;
        ++nStudy;

        //Store particles momentum in a file
        file << pythia.event[i].px() << " "<< setw(11) << pythia.event[i].py() 
             << " " << setw(11) << pythia.event[i].pz()
             << " " << setw(11) << pythia.event[i].e() << " " << setw(11) <<  endl;

      // Store momenta. Use energy component for absolute momentum.
      Vec4 pNow = pythia.event[i].p();
      pNow.e(pNow.pAbs());
      pSum += pNow;
      pOrder.push_back(pNow);
      }

  // Very low multiplicities (0 or 1) not considered.
    if (nStudy < 2) {
      if (nFew < 1) cout << " PYTHIA Error in "
        << "Thrust::analyze: too few particles" << endl;
      ++nFew;
      return false;
    }

    aRef = lin.eventAxis(1);
    //aRef /=aRef.pAbs();
    double dSum = 0.0;
  // Add all momenta with sign; two choices for each reference particle.
    for (int i = 0; i < nStudy; ++i){
      if (dot3(pOrder[i], aRef) > 0.){
        dSum += dot3(pOrder[i], aRef); 
        pPart += pOrder[i];
      }
      //else                            pPart -= pOrder[i];
    }

    //Thrust axis and value
    eVal1 = dSum / pSum.e();
    firstVal = eVal1;
    //eVal1 = pPart.e() / pSum.e();
    eVec1 = pPart / pPart.e();
    eVec1.e(0.);

    //cout << "eVal1: " << 2*eVal1  << " tMax: "<< 2*tMax << endl;

    while(tMax <= eVal1){
      //std::cout << "-----------Entering new iteration!---------" << std::endl;

      double newSum = 0.0;
      Vec4 newPart = 0.;
      pPart /= pPart.pAbs();

       // Add all momenta with sign; two choices for each reference particle.
      for (int i = 0; i < nStudy; ++i){
        if (dot3(pOrder[i], pPart) > 0.){
        newSum += dot3(pOrder[i], pPart); 
        newPart += pOrder[i];
        }
      }
      pPart = 0.;
      pPart = newPart;
      tMax = newSum / pSum.e();

      //std:: cout << "New value: " << 2*tMax << std::endl;

      if(tMax > eVal1){
        nIter++;
        eVal1 = tMax;
      //  cout << "Next iteration" << endl;
      }
      else break; 
 
    } // end of while loop

    //double diff = 0.;
    //diff = 2*tMax - thr.thrust();  
    //if ( diff < 1e-3)
    //  diff = 0.0;

    //Store particles momentum in a file
    filer << nStudy << " " << 2*tMax << " "<< setw(11) << thr.thrust() 
             << " " << setw(11) << nIter << " " << setw(11) << 2*firstVal << " "<< endl;     

    cout << "----------Alternative thrust value ------------\n"
         << "First value: " << setw(11) << 2*firstVal << setw(11)
         << "Value " << setw(11) << 2*tMax << setw(11) 
         << "Thrust: " << setw(11) << thr.thrust() 
         << "nIter: "<< setw(11) << nIter<< endl; 

  // End of event loop. Statistics. Output histograms.
  }
  pythia.stat();
  
  // Done.
  return 0;
}
