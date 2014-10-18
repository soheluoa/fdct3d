
/*  =================================================================================
 *  Name        : Fista.cpp
 *  Author      : Sohel Bhuiyan
 *  Version     : 1.0
 *  Purpose     : To reconstruct 3D/2D seismic data
 *  Date        : September 12, 2014
 *  Affiliation : University of Alberta, Physics department, (SAIG)
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
   fsp.maxItr = 10;
   fsp.eigenItr = 10;
   fsp.stat = 0;
   fsp.lambda = 0.0;
   fsp.inDataFileName = "/home/entropy/workspace/fdct3d/src/Data/data_gain.bin";
   fsp.sampleMatFileName = "/home/entropy/workspace/fdct3d/src/Data/sample_mat.bin";
   fsp.outDataFileName = "/home/entropy/workspace/fdct3d/src/Data/reconData.bin";
   fsp.testDataFileName = "/home/entropy/workspace/fdct3d/src/Data/testData.bin";

   /* Data initialisation for 3D FFT ****************************************************/
   params.ac = 0;

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

   CpxNumTns init(fsp.n3,fsp.n2,fsp.n1);

   for(int i=0; i<fsp.n3; i++)
 	  for(int j=0; j<fsp.n2; j++)
 		  for(int k=0; k<fsp.n1; k++)
 			  init(i,j,k) = cpx(1.0, 0.0);

   // Initialise the parameters for 3D curvelet transformation ***************************
   fdct3d_param(fsp.n3,fsp.n2,fsp.n1,params.nbscales,params.nbdstz_coarse,params.ac,params.fxs,params.fys,params.fzs, params.nxs,params.nys,params.nzs);

   fsc.readSampleMat(fsp.sampleMatFileName, fsp.n1, fsp.n2, fsp.n3, params.samplMat);

   fsp.alpha = 1; //*fsc.powerEigen(init, &fsp,  &params);

   fsp.thresh = fsp.lambda/(2*fsp.alpha);

   params.cellStruct.clear();

   //* Executing FISTA optimisation code to reconstruct seismic data *********************/
   start = clock();
   fsc.Reconstruct(&fsp, &params, &fsc);
   finish = clock();

   std::cout << "Time: " << (finish-start)/double(CLOCKS_PER_SEC) << " Seconds " <<std::endl;

   /** Clearing the memory space **********************************************************/
   std::cout << "Cleared the memory space." <<std::endl;
   fsp.reset();
   params.reset();
   return 0;
}

