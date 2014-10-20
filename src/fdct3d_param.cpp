/*
    FDCT3D (Fast 3d Curvelet Transform)
	Written by Lexing Ying

	Updated by: Sohel Bhuiyan (2014)
				University of Alberta
	   	   	   	Upgraded from FFTW 2.1.5 to FFTW 3.3.3
	   	   	    Solved the memory leaking problem in FFTW 3.3.3
	   	   	    Plugged in with FISTA Algorithm
*/

#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

//----------------------------------------
int fdct3d_param_center(float L1, float L2, float L3, int s, vector< vector<float> >& fxs, vector< vector<float> >& fys, vector< vector<float> >& fzs, vector< vector<int   > >& nxs, vector< vector<int   > >& nys, vector< vector<int   > >& nzs)
{
  vector<float>& fx = fxs[s];  vector<float>& fy = fys[s];  vector<float>& fz = fzs[s];
  vector<int>& nx = nxs[s];  vector<int>& ny = nys[s];  vector<int>& nz = nzs[s];
  
  int S1, S2, S3;	 int F1, F2, F3;	 float R1, R2, R3;
  fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
  fx.resize(1);  fy.resize(1);  fz.resize(1);
  nx.resize(1);  ny.resize(1);  nz.resize(1);
  
  fx[0] = 0;  fy[0] = 0;  fz[0] = 0;
  nx[0] = S1;  ny[0] = S2;  nz[0] = S3;
  
  return 0;
}

int fdct3d_param_wavelet(float L1, float L2, float L3, int s, int N1,int N2,int N3,
							 vector< vector<float> >& fxs, vector< vector<float> >& fys, vector< vector<float> >& fzs,
							 vector< vector<int   > >& nxs, vector< vector<int   > >& nys, vector< vector<int   > >& nzs)
{
  vector<float>& fx = fxs[s];  vector<float>& fy = fys[s];  vector<float>& fz = fzs[s];
  vector<int>& nx = nxs[s];  vector<int>& ny = nys[s];  vector<int>& nz = nzs[s];
  
  fx.resize(1);  fy.resize(1);  fz.resize(1);
  nx.resize(1);  ny.resize(1);  nz.resize(1);
  
  fx[0] = 0;  fy[0] = 0;  fz[0] = 0;
  nx[0] = N1;  ny[0] = N2;  nz[0] = N3;

  return 0;
}

//----------------------------------------
int fdct3d_param_angles(float L1, float L2, float L3, int s, int nd,
							vector< vector<float> >& fxs, vector< vector<float> >& fys, vector< vector<float> >& fzs,
							vector< vector<int   > >& nxs, vector< vector<int   > >& nys, vector< vector<int   > >& nzs)
{
  vector<float>& fx = fxs[s];  vector<float>& fy = fys[s];  vector<float>& fz = fzs[s];
  vector<int>& nx = nxs[s];  vector<int>& ny = nys[s];  vector<int>& nz = nzs[s];
  
  int nbw = 6 * nd * nd;
  fx.resize(nbw);  fy.resize(nbw);  fz.resize(nbw);
  nx.resize(nbw);  ny.resize(nbw);  nz.resize(nbw);
  
  int nf = 6;
  int wcnt = 0;
  int S1, S2, S3;	 int F1, F2, F3;	 float R1, R2, R3;	 fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
  float W1 = L1/nd;  float W2 = L2/nd;  float W3 = L3/nd;

  //face 0: x,y,z
  for(int h=0; h<nd; h++) { //(y first z second)
	 for(int g=0; g<nd; g++) {
		float xs = R1/4-(W1/2)/4;		float xe = R1;
		float ys = -R2 + (2*g-1)*W2/2;		float ye = -R2 + (2*g+3)*W2/2;
		float zs = -R3 + (2*h-1)*W3/2;		float ze = -R3 + (2*h+3)*W3/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		float tx, ty, tz;
		tx = R1/2;		  ty = (-R2 + (g+0.5)*W2)/2;		  tz = (-R3 + (h+0.5)*W3)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }
  //face 1: y z x
  for(int f=0; f<nd; f++) {
	 for(int h=0; h<nd; h++) {
		float ys = R2/4-(W2/2)/4;		  float ye = R2;
		float zs = -R3 + (2*h-1)*W3/2;		  float ze = -R3 + (2*h+3)*W3/2;
		float xs = -R1 + (2*f-1)*W1/2;		  float xe = -R1 + (2*f+3)*W1/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		float tx, ty, tz;
		ty = R2/2;		  tz = (-R3 + (h+0.5)*W3)/2;		  tx = (-R1 + (f+0.5)*W1)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }
  //face 2: z,x,y
  for(int g=0; g<nd; g++) {
	 for(int f=0; f<nd; f++) {
		float zs = R3/4-(W3/2)/4;		float ze = R3;
		float xs = -R1 + (2*f-1)*W1/2;		float xe = -R1 + (2*f+3)*W1/2;
		float ys = -R2 + (2*g-1)*W2/2;		float ye = -R2 + (2*g+3)*W2/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		float tx, ty, tz;
		tz = R3/2;		  tx = (-R1 + (f+0.5)*W1)/2;		  ty = (-R2 + (g+0.5)*W2)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }
  //face 3: -x,-y,-z
  for(int h=nd-1; h>=0; h--) {
	 for(int g=nd-1; g>=0; g--) {
		float xs = -R1;		  float xe = -R1/4+(W1/2)/4;
		float ys = -R2 + (2*g-1)*W2/2;		float ye = -R2 + (2*g+3)*W2/2;
		float zs = -R3 + (2*h-1)*W3/2;		float ze = -R3 + (2*h+3)*W3/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		float tx, ty, tz;
		tx = -R1/2;		  ty = (-R2 + (g+0.5)*W2)/2;		  tz = (-R3 + (h+0.5)*W3)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }

  //face 4: -y,-z,-x
  for(int f=nd-1; f>=0; f--) {
	 for(int h=nd-1; h>=0; h--) {
		float ys = -R2;		  float ye = -R2/4+(W2/2)/4;
		float zs = -R3 + (2*h-1)*W3/2;		  float ze = -R3 + (2*h+3)*W3/2;
		float xs = -R1 + (2*f-1)*W1/2;		  float xe = -R1 + (2*f+3)*W1/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		float tx, ty, tz;
		ty = -R2/2;		  tz = (-R3 + (h+0.5)*W3)/2;		  tx = (-R1 + (f+0.5)*W1)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }
  
  //face 5: -z,-x,-y
  for(int g=nd-1; g>=0; g--) {
	 for(int f=nd-1; f>=0; f--) {
		float zs = -R3;		  float ze = -R3/4+(W3/2)/4;
		float xs = -R1 + (2*f-1)*W1/2;		float xe = -R1 + (2*f+3)*W1/2;
		float ys = -R2 + (2*g-1)*W2/2;		float ye = -R2 + (2*g+3)*W2/2;
		int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
		float tx, ty, tz;
		tz = -R3/2;		  tx = (-R1 + (f+0.5)*W1)/2;		  ty = (-R2 + (g+0.5)*W2)/2;
		fx[wcnt] = tx;		  fy[wcnt] = ty;		  fz[wcnt] = tz;
		nx[wcnt] = xn;		  ny[wcnt] = yn;		  nz[wcnt] = zn;
		wcnt++;
	 }
  }
  
  return 0;
}

//----------------------------------------
int fdct3d_param(int N1, int N2, int N3, int nbscales, int nbdstz_coarse, int ac,
				  vector< vector<float> >& fxs, vector< vector<float> >& fys, vector< vector<float> >& fzs,
				  vector< vector<int   > >& nxs, vector< vector<int   > >& nys, vector< vector<int   > >& nzs)
{
  fxs.resize(nbscales);  fys.resize(nbscales);  fzs.resize(nbscales);
  nxs.resize(nbscales);  nys.resize(nbscales);  nzs.resize(nbscales);
  
  int L = nbscales;

  if(ac==1) {
	 {
		int s = 0;
		float L1 = 4.0*N1/3.0 / pow2(L-1-s);	 float L2 = 4.0*N2/3.0 / pow2(L-1-s);	 float L3 = 4.0*N3/3.0 / pow2(L-1-s);
		fdct3d_param_center(L1,L2,L3,s, fxs,fys,fzs, nxs,nys,nzs);
	 }
	 for(int s=1; s<L; s++) {
		float L1 = 4.0*N1/3.0 / pow2(L-1-s);	 float L2 = 4.0*N2/3.0 / pow2(L-1-s);	 float L3 = 4.0*N3/3.0 / pow2(L-1-s);
		int nd = nbdstz_coarse * pow2(s/2);
		fdct3d_param_angles(L1,L2,L3,s, nd, fxs,fys,fzs, nxs,nys,nzs);
	 }
  } else {
	 {
		int s = 0;
		float L1 = 4.0*N1/3.0 / pow2(L-1-s);	 float L2 = 4.0*N2/3.0 / pow2(L-1-s);	 float L3 = 4.0*N3/3.0 / pow2(L-1-s);
		fdct3d_param_center(L1,L2,L3,s, fxs,fys,fzs, nxs,nys,nzs);
	 }
	 for(int s=1; s<L-1; s++) {
		float L1 = 4.0*N1/3.0 / pow2(L-1-s);	 float L2 = 4.0*N2/3.0 / pow2(L-1-s);	 float L3 = 4.0*N3/3.0 / pow2(L-1-s);
		int nd = nbdstz_coarse * pow2(s/2);

		// Define angles for 6 different faces of a cube ***************************************
		fdct3d_param_angles(L1,L2,L3,s, nd, fxs,fys,fzs, nxs,nys,nzs);
	 }
	 {
		int s = L-1;
		float L1 = 4.0*N1/3.0 / pow2(L-1-s);	 float L2 = 4.0*N2/3.0 / pow2(L-1-s);	 float L3 = 4.0*N3/3.0 / pow2(L-1-s);
		fdct3d_param_wavelet(L1,L2,L3,s, N1,N2,N3, fxs,fys,fzs, nxs,nys,nzs);
	 }
  }
  
  return 0;
}
