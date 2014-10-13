/*
    FDCT3D (Fast 3d Curvelet Transform)
	Written by Lexing Ying

	Updated by: Sohel Bhuiyan (2014)
			    University of Alberta
	            Upgrade from FFTW 2.1.5 to FFTW 3.3.3
	            Solved the memory leaking problem in FFTW 3.3.3
*/

#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

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
  time_t tm0, tm1;

  assert(argc==3);
  map<string, string> opts;
  optionsCreate(argv[2], opts);
  map<string,string>::iterator mi;
  
  int m;
  mi = opts.find("-m"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>m; }
  int n;
  mi = opts.find("-n"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>n; }
  int p;
  mi = opts.find("-p"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>p; }
  
  int nbscales;
  mi = opts.find("-nbscales"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>nbscales; }
  
  int nbdstz_coarse;
  mi = opts.find("-nbdstz_coarse"); assert(mi!=opts.end());
  { istringstream ss((*mi).second); ss>>nbdstz_coarse; }

  int ac = 0;

  std::cout << m << " " << n << " " << p << std::endl;

  srand48( (long)time(NULL) );

  // Class of NumTns *********************************************************************************
  CpxNumTns x(m,n,p);

  for(int i=0; i<m; i++)
	for(int j=0; j<n; j++)
		for(int k=0; k<p; k++)
		{
			x(i,j,k) =  cpx(drand48()-0.5,drand48());
			std::cout << x(i,j,k)  <<std::endl;
		}

  tm0 = time(NULL);

  vector< vector<double> > fxs,fys,fzs;
  vector< vector<int> > nxs,nys,nzs;
  fdct3d_param(m, n, p, nbscales, nbdstz_coarse,ac,fxs,fys,fzs,nxs,nys,nzs);
  tm1 = time(NULL);
  cout<<"fdct3d_param "<<difftime(tm1,tm0)<<" seconds"<<endl;  tm0 = tm1;
  
  // Vector of class NumTns with complex double type ***************************************************
  vector< vector<CpxNumTns> > c;

  fdct3d_forward(m, n, p, nbscales, nbdstz_coarse, ac, x, c);
  tm1 = time(NULL);
  cout<<"fdct3d_forward "<<difftime(tm1,tm0)<<" seconds"<<endl;  tm0 = tm1;
  
  CpxNumTns newx(x);
  clear(newx);
  fdct3d_inverse(m, n, p, nbscales, nbdstz_coarse, ac, c, newx);

  //std::cout << (c(0,0)) << std::endl;
  tm1 = time(NULL);
  cout<<"fdct3d_inverse "<<difftime(tm1,tm0)<<" seconds"<<endl;
  tm0 = tm1;  //cerr<<energy(newx)<<endl;
  
  double mv = 0.0;
  double sum = 0.0;
  for(int i=0; i<m; i++)
 	for(int j=0; j<n; j++)
		for(int k=0; k<p; k++)
		{
			mv = max(mv, abs(newx(i,j,k)-x(i,j,k)));
		    sum += abs(newx(i,j,k)-x(i,j,k));
		}

  cerr<<"max error "<<mv<<endl;
  cerr<<"Total error "<<sum<<endl;

  return 0;
}
