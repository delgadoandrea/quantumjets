// Code generated to use input particles clustered by qubo formulation as seed axis to calculate thrust
#include "Pythia8/Pythia.h"
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;
using namespace Pythia8;

const int SIZE = 100;
double px[SIZE];
double py[SIZE];
double pz[SIZE];
double e[SIZE];

vector<int> solution(string fname){
	int num = 0;
	vector<int> solution;

	string inFileName = fname;
	ifstream inFile;
	inFile.open(inFileName.c_str());

	if (inFile.is_open()){
		while(!inFile.eof()){

			int value = 0;
			inFile >> value;
			solution.push_back(value);
			//std::cout << value << std::endl;
			++num;
		}
		inFile.close();
	}
	else{
		cerr << "Can't find file with D-Wave solutions: " << inFileName << endl;
	}

	return solution;

}

double nref(string fname, string fsol){
	int num = 0;
	int nIter = 0;
	Vec4 pNow, aRef, pSum;
	vector<Vec4> pOrder;

	string inFileName = fname;//"/home/andrea/pythiaEvents/Event_0.dat";
	string inFileNameSolutions = fsol;


	ifstream inFile;
	inFile.open(inFileName.c_str());

	if (inFile.is_open()){
		while(!inFile.eof()){
			
			inFile >> px[num];
			inFile >> py[num];
			inFile >> pz[num];
			inFile >> e[num];

			pNow.px(px[num]);
			pNow.py(py[num]);
			pNow.pz(pz[num]);
			pNow.e(e[num]);

			//Use energy component for absolute momentum
			pNow.e(pNow.pAbs());

			pSum += pNow;
			pOrder.push_back(pNow);

			//std::cout << px[num] << " " << py[num] << " " << pz[num] << " " << e[num] << std::endl;

			++num;
		}
		inFile.close();
	}
	else{
		cerr << "Can't find input file: " << inFileName << endl;
	}

	vector<int> solVec;
	solVec = solution(fsol);
	
	for(int i = 0; i < num; i++) if(solVec[i] == 1) aRef += pOrder[i];

	//cout << "Sum of partition: " << aRef.pAbs() << endl;
	//cout << "Total sum e: " << pSum.e() << endl;
	//cout << "D-Wave Thrust: " << 2*aRef.pAbs()/pSum.e() << setprecision(8) << endl;

	//-------------------Start of thrust calculation using seed axis --------
    double eVal1 = 0.;
    double tMax = 0.0;
    double dSum = 0.0;
	Vec4 pPart, eVec1;
  	// Add all momenta with sign; two choices for each reference particle.
    for (int i = 0; i < num; ++i){
     	if (dot3(pOrder[i], aRef) > 0.){
        	dSum += dot3(pOrder[i], aRef); 
        	pPart += pOrder[i];
      	}
    }
    //tMax = aRef.pAbs()/pSum.e();
    tMax = 0.0;

    //Thrust axis and value
    eVal1 = pPart.pAbs() / pSum.e();
    //cout << "------Seeded thrust: " << 2*eVal1 << setprecision(10)<< std::endl;
    eVec1 = pPart / pPart.e();
    eVec1.e(0.);

    //cout << "eVal1: " << 2*eVal1 << std::endl;//2*sqrt(eVal1)/pSum.e() << endl;

    while(tMax <= eVal1){
    	
    	//std::cout << "------Entered new iteration-------" <<std::endl;

    	double newSum = 0.0;
    	Vec4 newPart = 0.;
    	pPart /= pPart.pAbs();
    	//cout << pPart.e() << endl;

       // Add all momenta with sign; two choices for each reference particle.
    	for (int i = 0; i < num; ++i){
        	if (dot3(pOrder[i], pPart) > 0.){
        		newSum += dot3(pOrder[i], pPart); 
        		newPart += pOrder[i];
        	}
      	}
      
      	pPart.reset();
      	pPart = newPart;
      	//cout << pPart.e() << endl;
      	tMax = newSum / pSum.e();

      	//std:: cout << "------New value: " << 2*tMax << setprecision(10) << std::endl;

    	if(tMax > eVal1){
    		nIter++;
        	eVal1 = tMax;
        	//cout << "Next iteration" << endl;
        	
      	}
      	else break; 



    } // end of while loop

    //cout << "*****Number of iterations: " << nIter << endl;
    //double diff = 0.;
    //diff = 2*tMax - thr.thrust();  
    //if ( diff < 1e-3)
    //diff = 0.0;

            //Store particles momentum in a file
    //filer << diff << " "<< setw(11) << thr.thrust() 
    //         << " " << setw(11) <<  endl;    

    //filer << 2*tMax << " "<< setw(11) << nIter 
    //         << " " << setw(11) << endl;    
    if(nIter > 2){
    	cout << "----------Alternative thrust value ------------\n"
         	<< "Value " << setw(11) << 2*tMax << setw(11)  << " " << nIter 
         	<< " Event number: " << fname << endl; 
    }
    //cout << nIter << endl;
    return 2*tMax;

}

int main(){

	std::string path = "/home/andrea/pythiaEvents/";
	std::string fres = path + "DWaveSeedResults.dat";
	ofstream filer(fres);

	for (int iEvent = 0; iEvent < 1000; iEvent++){
		std::string fev    = path + "Event_" + std::to_string(iEvent) + ".dat";
		std::string fevsol = path + "Event_" + std::to_string(iEvent) + "_qsolutions.dat";

		double solution = nref(fev, fevsol);

		filer << solution << " "<< setw(11) << endl; 
	}
	return 0;
}


