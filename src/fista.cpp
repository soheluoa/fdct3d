/*
 * spgl1.cpp
 *
 *  Created on: Sep 12, 2014
 *      Author: entropy
 */

/*  =================================================================================
 *  Name        : Migration.cpp (Phase Shift Migration Plus Interpolation and others)
 *  Author      : Sohel Bhuiyan
 *  Version     : 0.1
 *  Purpose     : Generate the migrated image of the velocity model using PSPI method
 *  Date        : Feb 18, 2014 - implementation for PSPI poststack method
 *  Affiliation : University of Alberta, Physics department (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */

#include "fistaParams.hpp"
#include "fistaCore.hpp"
#include "PARAMS.hpp"
#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

#include <cmath>


int optionsCreate(const char* optfile, map<string,string>& options)
{

  options.clear();
  ifstream fin(optfile);
  fin.open("/home/entropy/workspace/fdct3d/src/params.txt");
  std::cout << optfile << std::endl;
  std::cout <<fin.fail()  << std::endl;
  assert(fin.good());

  string name;
  fin>>name;

  while(fin.good()) {
	 char cont[100];
	 fin.getline(cont, 99);
	 options[name] = string(cont);
	 fin>>name;
  }
  fin.close();
  return 0;
}

int main(int argc, char** argv)
{
   clock_t start, finish;

   /** Creating the objects of the dependence classes ************************************/
   fistaParams fsp;
   fistaCore fsc;
   PARAMS params;

   /* Data initialisation */
   fsp.bpTol = 0.000001;
   fsp.maxItr = 200;
   fsp.stat = 0;
   fsp.n1 = 500;
   fsp.n2 = 36;
   fsp.n3 =	100;
   fsp.lambda = 1.05;
   fsp.alpha = 2.02;
   fsp.inDataFileName = "/home/entropy/workspace/fista/src/data_gain.bin";
   fsp.sampleMatFileName = "/home/entropy/workspace/fista/src/sample_mat.bin";
   fsp.outDataFileName = "reconData.bin";

   /* Data initialisation for 3D FFT ****************************************************/
   params.Nw = pow(2,ceil(log10(fsp.n1)/log10(2)));
   params.Kx = pow(2,ceil(log10(fsp.n2)/log10(2)));
   params.Ky = pow(2,ceil(log10(fsp.n3)/log10(2)));
   params.ac = 0;

   fsp.thresh = fsp.lambda/(2*fsp.alpha);

   assert(argc==3);
   map<string, string> opts;
   optionsCreate(argv[2], opts);
   map<string,string>::iterator mi;

   mi = opts.find("-m");
   assert(mi!=opts.end());
   {istringstream ss((*mi).second); ss>>fsp.n1;}

   mi = opts.find("-n"); assert(mi!=opts.end());
   {istringstream ss((*mi).second); ss>>fsp.n2;}

   mi = opts.find("-p"); assert(mi!=opts.end());
   {istringstream ss((*mi).second); ss>>fsp.n3;}

   mi = opts.find("-nbscales"); assert(mi!=opts.end());
   {istringstream ss((*mi).second); ss>>params.nbscales;}

   mi = opts.find("-nbdstz_coarse"); assert(mi!=opts.end());
   {istringstream ss((*mi).second); ss>>params.nbdstz_coarse;}

   /** Executing FISTA optimisation code to reconstruct seismic data *********************/
   start = clock();
   fsc.Execute(&fsp, &params, &fsc);
   finish = clock();

   std::cout << "Time: " << (finish-start)/double(CLOCKS_PER_SEC) << " Seconds " <<std::endl;

   /** Clearing the memory space **********************************************************/
   std::cout << "Cleared the memory space." <<std::endl;
   fsp.reset();
   params.reset();
   return 0;
}
