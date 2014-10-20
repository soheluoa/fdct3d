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

//-----------------------------------------------------------------------------------------------------------------------------
int fdct3d_forward_angles(float L1, float L2, float L3, int s, int nd, CpxOffTns& O, vector< vector<CpxNumTns> >& C, std::vector<std::vector<std::vector<int> > > &cellStruct)
{
  //allocate space
  vector<CpxNumTns>& csc = C[s];
  cellStruct.push_back(vector<vector<int> > ());

  csc.resize(6*nd*nd);
  
  int nf = 6;
  int wcnt = 0;
  int S1, S2, S3;
  int F1, F2, F3;
  float R1, R2, R3;
  fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);
  DblOffVec big1(S1);  fdct3d_lowpass(L1, big1);
  DblOffVec big2(S2);  fdct3d_lowpass(L2, big2);
  DblOffVec big3(S3);  fdct3d_lowpass(L3, big3);
  
  float Lh1 = L1/2;  float Lh2 = L2/2;  float Lh3 = L3/2;
  int Sh1, Sh2, Sh3;   int Fh1, Fh2, Fh3;	 float Rh1, Rh2, Rh3;
  fdct3d_rangecompute(Lh1, Lh2, Lh3, Sh1, Sh2, Sh3, Fh1, Fh2, Fh3, Rh1, Rh2, Rh3);
  DblOffVec sma1(S1);  fdct3d_lowpass(Lh1, sma1);
  DblOffVec sma2(S2);  fdct3d_lowpass(Lh2, sma2);
  DblOffVec sma3(S3);  fdct3d_lowpass(Lh3, sma3);
  
  float W1 = L1/nd;  float W2 = L2/nd;  float W3 = L3/nd;
  
  typedef pair<int,int> intpair;
  typedef pair<int, intpair> inttriple;
  map<inttriple, fftwf_plan> planmap;

  //face 0: x,y,z
  for(int h=0; h<nd; h++) { //(y first z second)
	 for(int g=0; g<nd; g++) {
	   float xs = R1/4-(W1/2)/4;
	   float xe = R1;
	   float ys = -R2 + (2*g-1)*W2/2;
	   float ye = -R2 + (2*g+3)*W2/2;
	   float zs = -R3 + (2*h-1)*W3/2;
	   float ze = -R3 + (2*h+3)*W3/2;
	   int xn = int(ceil(xe-xs));
	   int yn = int(ceil(ye-ys));
	   int zn = int(ceil(ze-zs));
	   float thts, thtm, thte; //y to x

	   if(g==0) {
		 thts = atan2(-1.0, 1.0-1.0/nd);
		 thtm = atan2(-1.0+1.0/nd, 1.0);
		 thte = atan2(-1.0+3.0/nd, 1.0);
	   } else if(g==nd-1) {
		 thts = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);
		 thtm = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);
		 thte = atan2(1.0, 1.0-1.0/nd);
	   } else {
		 thts = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);
		 thtm = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);
		 thte = atan2(-1.0+(2.0*g+3.0)/nd, 1.0);
	   }
	   float phis, phim, phie; //z to x

	   if(h==0) {
		 phis = atan2(-1.0, 1.0-1.0/nd);
		 phim = atan2(-1.0+1.0/nd, 1.0);
		 phie = atan2(-1.0+3.0/nd, 1.0);
	   } else if(h==nd-1) {
		 phis = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);
		 phim = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);
		 phie = atan2(1.0, 1.0-1.0/nd);
	   } else {
		 phis = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);
		 phim = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);
		 phie = atan2(-1.0+(2.0*h+3.0)/nd, 1.0);
	   }

	   int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2; //half
	   float R21 = R2/R1;		  float R31 = R3/R1;
	   CpxOffTns wpdata(xn,yn,zn);

	   for(int xcur=(int)ceil(xs); xcur<xe; xcur++)
	   {
		 int yfm = (int)ceil( max(-R2, R21*xcur*tan(thts)) );
		 int yto = (int)floor( min(R2, R21*xcur*tan(thte)) );
		 int zfm = (int)ceil( max(-R3, R31*xcur*tan(phis)) );
		 int zto = (int)floor( min(R3, R31*xcur*tan(phie)) );

		 for(int ycur=yfm; ycur<=yto; ycur++)
		   for(int zcur=zfm; zcur<=zto; zcur++)
		   {
			 int tmpx = xcur%xn;
			 if(tmpx<-xh)
				 tmpx+=xn;
			 if(tmpx>=-xh+xn)
				 tmpx-=xn;
			 int tmpy = ycur%yn;
			 if(tmpy<-yh)
				 tmpy+=yn;
			 if(tmpy>=-yh+yn)
				 tmpy-=yn;
			 int tmpz = zcur%zn;
			 if(tmpz<-zh)
				 tmpz+=zn;
			 if(tmpz>=-zh+zn)
				 tmpz-=zn;
			 float ss = sma1(xcur)*sma2(ycur)*sma3(zcur);
			 float bb = big1(xcur)*big2(ycur)*big3(zcur);

			 wpdata(tmpx, tmpy, tmpz) = O(xcur,ycur,zcur) * bb * sqrtf(1.0-ss*ss);

			 float thtcur = atan2(ycur/R2, xcur/R1);
			 float phicur = atan2(zcur/R3, xcur/R1);

			 float glbpou;
			 fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou);

			 float wtht;
			 if(thtcur<thtm)
			 {
			   if(g==0)
				   wtht = 1;
			   else
			   {
				   float l,r;
				   fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);
				   wtht = l;
			   }
			 }
			 else
			 {
			   if(g==nd-1)
				   	wtht = 1;
			   else
			   {
				   float l,r;
				   fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);
				   wtht = r;
			   }
			 }

			 float wphi;
			 if(phicur<phim) {
			   if(h==0)						wphi = 1;
			   else
			   {
				   float l,r;
				   fdct3d_window( (phicur-phis)/(phim-phis), l, r);
				   wphi = l;
			   }
			 } else {
			   if(h==nd-1)
				   wphi = 1;
			   else
			   {
				   float l,r;
				   fdct3d_window( (phicur-phim)/(phie-phim), l, r);
				   wphi = r;
			   }
			 }

			 float pou = glbpou * wtht * wphi;
			 wpdata(tmpx, tmpy, tmpz) *= pou;
		   }
	   } //xcur

	   CpxNumTns tpdata(xn,yn,zn);

	   fdct3d_ifftshift(xn,yn,zn,wpdata,tpdata);

	   fftwf_plan p = NULL;
	   p = fftwf_plan_dft_3d(zn, yn, xn, (fftwf_complex*)tpdata.data(), (fftwf_complex*)tpdata.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
	   fftwf_execute(p);		  		//cerr<<"wedge s"<<endl;

	   float sqrtprod = sqrt(float(xn*yn*zn));
	   for(int i=0; i<xn; i++)
		   for(int j=0; j<yn; j++)
			   for(int k=0; k<zn; k++)
				   tpdata(i,j,k) /= sqrtprod;

	   fftwf_destroy_plan(p);
	   cellStruct[s].push_back(vector<int> ());
	   cellStruct[s][wcnt].push_back(xn);
	   cellStruct[s][wcnt].push_back(yn);
	   cellStruct[s][wcnt].push_back(zn);

	   csc[wcnt] = tpdata;
	   wcnt++;
	 }
  } //end of face

  //face 1. y z x
  for(int f=0; f<nd; f++) {
	 for(int h=0; h<nd; h++) {
	   
	   float ys = R2/4-(W2/2)/4;		      float ye = R2;
	   float zs = -R3 + (2*h-1)*W3/2;		  float ze = -R3 + (2*h+3)*W3/2;
	   float xs = -R1 + (2*f-1)*W1/2;		  float xe = -R1 + (2*f+3)*W1/2;
	   int xn = int(ceil(xe-xs));		      int yn = int(ceil(ye-ys));
	   int zn = int(ceil(ze-zs));
	   float thts, thtm, thte; //z to y

	   if(h==0) {
		 thts = atan2(-1.0, 1.0-1.0/nd);
		 thtm = atan2(-1.0+1.0/nd, 1.0);
		 thte = atan2(-1.0+3.0/nd, 1.0);
	   }
	   else if(h==nd-1)
	   {
		 thts = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);
		 thtm = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);
		 thte = atan2(1.0, 1.0-1.0/nd);
	   }
	   else
	   {
		 thts = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);
		 thtm = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);
		 thte = atan2(-1.0+(2.0*h+3.0)/nd, 1.0);
	   }
	   float phis, phim, phie; //z to x

	   if(f==0) {
		 phis = atan2(-1.0, 1.0-1.0/nd);
		 phim = atan2(-1.0+1.0/nd, 1.0);
		 phie = atan2(-1.0+3.0/nd, 1.0);
	   }
	   else if(f==nd-1) {
		 phis = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);
		 phim = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);
		 phie = atan2(1.0, 1.0-1.0/nd);
	   } else {
		 phis = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);
		 phim = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);
		 phie = atan2(-1.0+(2.0*f+3.0)/nd, 1.0);
	   }

	   int xh = xn/2;
	   int yh = yn/2;
	   int zh = zn/2;
	   float R32 = R3/R2;
	   float R12 = R1/R2;

	   CpxOffTns wpdata(xn,yn,zn);

	   for(int ycur=(int)ceil(ys); ycur<ye; ycur++)
	   {
		 int zfm = (int)ceil( max(-R3, R32*ycur*tan(thts)) );
		 int zto = (int)floor( min(R3, R32*ycur*tan(thte)) );
		 int xfm = (int)ceil( max(-R1, R12*ycur*tan(phis)) );
		 int xto = (int)floor( min(R1, R12*ycur*tan(phie)) );

		 for(int zcur=zfm; zcur<=zto; zcur++)
		   for(int xcur=xfm; xcur<=xto; xcur++)
		   {
			 int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
			 int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
			 int tmpz = zcur%zn;				  if(tmpz<-zh) tmpz+=zn;				  if(tmpz>=-zh+zn) tmpz-=zn;

			 float ss = sma1(xcur)*sma2(ycur)*sma3(zcur);
			 float bb = big1(xcur)*big2(ycur)*big3(zcur);

			 wpdata(tmpx, tmpy, tmpz) = O(xcur,ycur,zcur) * bb * sqrtf(1.0-ss*ss);
			 float thtcur = atan2(zcur/R3, ycur/R2);
			 float phicur = atan2(xcur/R1, ycur/R2);
			 float glbpou;
			 fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou); //CHECK
			 float wtht;

			 if(thtcur<thtm) {
			   if(h==0)
				   wtht = 1;
			   else
			   {
				   float l,r;
				   fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);
				   wtht = l;
			   }
			 } else
			 {
			   if(h==nd-1)
				   wtht = 1;
			   else {
				   float l,r;
				   fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);
				   wtht = r;
			   }
			 }
			 float wphi;
			 if(phicur<phim) {
			   if(f==0)						wphi = 1;
			   else {
				   float l,r;
				   fdct3d_window( (phicur-phis)/(phim-phis), l, r);
				   wphi = l;
			   }
			 } else {
			   if(f==nd-1)
				   wphi = 1;
			   else {
				   float l,r;
				   fdct3d_window( (phicur-phim)/(phie-phim), l, r);
				   wphi = r;
			   }
			 }
			 float pou = glbpou * wtht * wphi;
			 wpdata(tmpx, tmpy, tmpz) *= pou;
		   }
	   } //ycur
	   CpxNumTns tpdata(xn,yn,zn);

	   fdct3d_ifftshift(xn,yn,zn,wpdata,tpdata);

	   fftwf_plan p = NULL;
	   p = fftwf_plan_dft_3d(zn, yn, xn, (fftwf_complex*)tpdata.data(), (fftwf_complex*)tpdata.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
	   fftwf_execute(p);		  //cerr<<"wedge s"<<endl;

	   float sqrtprod = sqrt(float(xn*yn*zn));
	   for(int i=0; i<xn; i++)
		   for(int j=0; j<yn; j++)
			   for(int k=0; k<zn; k++)
				   tpdata(i,j,k) /= sqrtprod;

	   fftwf_destroy_plan(p);

	   cellStruct[s].push_back(vector<int> ());
	   cellStruct[s][wcnt].push_back(xn);
	   cellStruct[s][wcnt].push_back(yn);
	   cellStruct[s][wcnt].push_back(zn);

	   csc[wcnt] = tpdata;
	   wcnt++;
	 }
  }//end of face

  //face 2. z x y
  for(int g=0; g<nd; g++) {
	 for(int f=0; f<nd; f++)
	 {
	   float zs = R3/4-(W3/2)/4;		float ze = R3;
	   float xs = -R1 + (2*f-1)*W1/2;		float xe = -R1 + (2*f+3)*W1/2;
	   float ys = -R2 + (2*g-1)*W2/2;		float ye = -R2 + (2*g+3)*W2/2;
	   int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
	   float thts, thtm, thte; //y to x
	   if(f==0) {
		 thts = atan2(-1.0, 1.0-1.0/nd);
		 thtm = atan2(-1.0+1.0/nd, 1.0);
		 thte = atan2(-1.0+3.0/nd, 1.0);
	   } else if(f==nd-1) {
		 thts = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);
		 thtm = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);
		 thte = atan2(1.0, 1.0-1.0/nd);
	   } else {
		 thts = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);
		 thtm = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);
		 thte = atan2(-1.0+(2.0*f+3.0)/nd, 1.0);
	   }
	   float phis, phim, phie; //z to x
	   if(g==0) {
		 phis = atan2(-1.0, 1.0-1.0/nd);
		 phim = atan2(-1.0+1.0/nd, 1.0);
		 phie = atan2(-1.0+3.0/nd, 1.0);
	   } else if(g==nd-1) {
		 phis = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);
		 phim = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);
		 phie = atan2(1.0, 1.0-1.0/nd);
	   } else {
		 phis = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);
		 phim = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);
		 phie = atan2(-1.0+(2.0*g+3.0)/nd, 1.0);
	   }
	   int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
	   float R13 = R1/R3;		  float R23 = R2/R3;//float R13 = float(F1)/float(F3);		  float R23 = float(F2)/float(F3);
	   CpxOffTns wpdata(xn,yn,zn);
	   for(int zcur=(int)ceil(zs); zcur<ze; zcur++) {
		 int xfm = (int)ceil( max(-R1, R13*zcur*tan(thts)) );
		 int xto = (int)floor( min(R1, R13*zcur*tan(thte)) );
		 int yfm = (int)ceil( max(-R2, R23*zcur*tan(phis)) );
		 int yto = (int)floor( min(R2, R23*zcur*tan(phie)) );

		 for(int xcur=xfm; xcur<=xto; xcur++)
		   for(int ycur=yfm; ycur<=yto; ycur++)
		   {
			 int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
			 int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
			 int tmpz = zcur%zn;				  if(tmpz<-zh) tmpz+=zn;				  if(tmpz>=-zh+zn) tmpz-=zn;
			 float ss = sma1(xcur)*sma2(ycur)*sma3(zcur);
			 float bb = big1(xcur)*big2(ycur)*big3(zcur);

			 wpdata(tmpx, tmpy, tmpz) = O(xcur,ycur,zcur) * bb * sqrtf(1.0-ss*ss);

			 float thtcur = atan2(xcur/R1, zcur/R3);
			 float phicur = atan2(ycur/R2, zcur/R3);

			 float glbpou;
			 fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou);

			 float wtht;
			 if(thtcur<thtm) {
			   if(f==0)						wtht = 1;
			   else {
				   float l,r;
				   fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);
				   wtht = l;
			   }
			 } else {
			   if(f==nd-1)
				   wtht = 1;
			   else {
				   float l,r;
				   fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);
				   wtht = r;
			   }
			 }
			 float wphi;
			 if(phicur<phim) {
			   if(g==0)						wphi = 1;
			   else {
				   float l,r;
				   fdct3d_window( (phicur-phis)/(phim-phis), l, r);
				   wphi = l;
			   }
			 } else {
			   if(g==nd-1)
				   wphi = 1;
			   else {
				   float l,r;
				   fdct3d_window( (phicur-phim)/(phie-phim), l, r);
				   wphi = r;
			   }
			 }
			 float pou = glbpou * wtht * wphi;
			 wpdata(tmpx, tmpy, tmpz) *= pou;
		   }
	   }//zcur
	   CpxNumTns tpdata(xn,yn,zn);
	   fdct3d_ifftshift(xn,yn,zn,wpdata,tpdata);

	   fftwf_plan p = NULL;
	   p = fftwf_plan_dft_3d(zn, yn, xn, (fftwf_complex*)tpdata.data(), (fftwf_complex*)tpdata.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
	   fftwf_execute(p);		  //cerr<<"wedge s"<<endl;

	   float sqrtprod = sqrt(float(xn*yn*zn));
	   for(int i=0; i<xn; i++)
		   for(int j=0; j<yn; j++)
			   for(int k=0; k<zn; k++)
				   tpdata(i,j,k) /= sqrtprod;

	   fftwf_destroy_plan(p);

	   cellStruct[s].push_back(vector<int> ());
	   cellStruct[s][wcnt].push_back(xn);
	   cellStruct[s][wcnt].push_back(yn);
	   cellStruct[s][wcnt].push_back(zn);

	   csc[wcnt] = tpdata;
	   wcnt++;
	 }
  }//end of face


///// ****************************************************
  //face 3: -x,-y,-z
  for(int h=nd-1; h>=0; h--) {
	 for(int g=nd-1; g>=0; g--) {
	   
	   float xs = -R1;		  float xe = -R1/4+(W1/2)/4;
	   float ys = -R2 + (2*g-1)*W2/2;		float ye = -R2 + (2*g+3)*W2/2;
	   float zs = -R3 + (2*h-1)*W3/2;		float ze = -R3 + (2*h+3)*W3/2;

	   int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));

	   float thts, thtm, thte; //y to x
	   if(g==0) {
		 thts = atan2(-1.0, 1.0-1.0/nd);
		 thtm = atan2(-1.0+1.0/nd, 1.0);
		 thte = atan2(-1.0+3.0/nd, 1.0);
	   }
	   else if(g==nd-1) {
		 thts = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);
		 thtm = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);
		 thte = atan2(1.0, 1.0-1.0/nd);
	   } else {
		 thts = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);
		 thtm = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);
		 thte = atan2(-1.0+(2.0*g+3.0)/nd, 1.0);
	   }
	   float phis, phim, phie; //z to x
	   if(h==0) {
		 phis = atan2(-1.0, 1.0-1.0/nd);
		 phim = atan2(-1.0+1.0/nd, 1.0);
		 phie = atan2(-1.0+3.0/nd, 1.0);
	   }
	   else if(h==nd-1)
	   {
		 phis = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);
		 phim = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);
		 phie = atan2(1.0, 1.0-1.0/nd);
	   } else
	   {
		 phis = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);
		 phim = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);
		 phie = atan2(-1.0+(2.0*h+3.0)/nd, 1.0);
	   }
	   int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
	   float R21 = R2/R1;		  float R31 = R3/R1;
	   CpxOffTns wpdata(xn,yn,zn);

	   for(int xcur=(int)ceil(xs); xcur<xe; xcur++) {
		 int yfm = (int)ceil( max(-R2, R21*(-xcur)*tan(thts)) );
		 int yto = (int)floor( min(R2, R21*(-xcur)*tan(thte)) );
		 int zfm = (int)ceil( max(-R3, R31*(-xcur)*tan(phis)) );
		 int zto = (int)floor( min(R3, R31*(-xcur)*tan(phie)) );

		 for(int ycur=yfm; ycur<=yto; ycur++)
		   for(int zcur=zfm; zcur<=zto; zcur++) {
			 int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
			 int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
			 int tmpz = zcur%zn;				  if(tmpz<-zh) tmpz+=zn;				  if(tmpz>=-zh+zn) tmpz-=zn;
			 float ss = sma1(xcur)*sma2(ycur)*sma3(zcur);
			 float bb = big1(xcur)*big2(ycur)*big3(zcur);

			 wpdata(tmpx, tmpy, tmpz) = O(xcur,ycur,zcur) * bb * sqrtf(1.0-ss*ss);
			 float thtcur = atan2(ycur/R2, (-xcur)/R1);
			 float phicur = atan2(zcur/R3, (-xcur)/R1);

			 float glbpou;
			 fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou);
			 float wtht;
			 if(thtcur<thtm) {
			   if(g==0)						wtht = 1;
			   else {
				   float l,r;
				   fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);
				   wtht = l;
			   }
			 } else {
			   if(g==nd-1)						wtht = 1;
			   else {
				   float l,r;
				   fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);
				   wtht = r;
			   }
			 }
			 float wphi;

			 if(phicur<phim) {
			   if(h==0)						wphi = 1;
			   else {
				   float l,r;
				   fdct3d_window( (phicur-phis)/(phim-phis), l, r);
				   wphi = l;
			   }
			 } else {
			   if(h==nd-1)						wphi = 1;
			   else {
				   float l,r;
				   fdct3d_window( (phicur-phim)/(phie-phim), l, r);
				   wphi = r;
			   }
			 }
			 float pou = glbpou * wtht * wphi;
			 wpdata(tmpx, tmpy, tmpz) *= pou;
		   }
	   } //xcur

	   CpxNumTns tpdata(xn,yn,zn);
	   fdct3d_ifftshift(xn,yn,zn,wpdata,tpdata);
	   fftwf_plan p = NULL;
	   p = fftwf_plan_dft_3d(zn, yn, xn, (fftwf_complex*)tpdata.data(), (fftwf_complex*)tpdata.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
	   fftwf_execute(p);		  //cerr<<"wedge s"<<endl;

	   float sqrtprod = sqrt(float(xn*yn*zn));
	   for(int i=0; i<xn; i++)
		   for(int j=0; j<yn; j++)
			   for(int k=0; k<zn; k++)
				   tpdata(i,j,k) /= sqrtprod;

	   fftwf_destroy_plan(p);

	   cellStruct[s].push_back(vector<int> ());
	   cellStruct[s][wcnt].push_back(xn);
	   cellStruct[s][wcnt].push_back(yn);
	   cellStruct[s][wcnt].push_back(zn);

	   csc[wcnt] = tpdata;
	   wcnt++;
	 }
  } //end of face
  //face 4. -y,-z,-x
  for(int f=nd-1; f>=0; f--) {
	for(int h=nd-1; h>=0; h--) {
	  
	  float ys = -R2;		  float ye = -R2/4+(W2/2)/4;
	  float zs = -R3 + (2*h-1)*W3/2;		  float ze = -R3 + (2*h+3)*W3/2;
	  float xs = -R1 + (2*f-1)*W1/2;		  float xe = -R1 + (2*f+3)*W1/2;
	  int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
	  float thts, thtm, thte; //z to y
	  if(h==0) {
		thts = atan2(-1.0, 1.0-1.0/nd);
		thtm = atan2(-1.0+1.0/nd, 1.0);
		thte = atan2(-1.0+3.0/nd, 1.0);
	  } else if(h==nd-1) {
		thts = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);
		thtm = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);
		thte = atan2(1.0, 1.0-1.0/nd);
	  } else {
		thts = atan2(-1.0+(2.0*h-1.0)/nd, 1.0);
		thtm = atan2(-1.0+(2.0*h+1.0)/nd, 1.0);
		thte = atan2(-1.0+(2.0*h+3.0)/nd, 1.0);
	  }
	  float phis, phim, phie; //z to x
	  if(f==0) {
		phis = atan2(-1.0, 1.0-1.0/nd);
		phim = atan2(-1.0+1.0/nd, 1.0);
		phie = atan2(-1.0+3.0/nd, 1.0);
	  } else if(f==nd-1) {
		phis = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);
		phim = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);
		phie = atan2(1.0, 1.0-1.0/nd);
	  } else {
		phis = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);
		phim = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);
		phie = atan2(-1.0+(2.0*f+3.0)/nd, 1.0);
	  }
	  int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
	  float R32 = R3/R2;		  float R12 = R1/R2;	  //float R32 = float(F3)/float(F2);		  float R12 = float(F1)/float(F2);

	  CpxOffTns wpdata(xn,yn,zn);

	  for(int ycur=(int)ceil(ys); ycur<ye; ycur++) {
		int zfm = (int)ceil( max(-R3, R32*(-ycur)*tan(thts)) );
		int zto = (int)floor( min(R3, R32*(-ycur)*tan(thte)) );
		int xfm = (int)ceil( max(-R1, R12*(-ycur)*tan(phis)) );
		int xto = (int)floor( min(R1, R12*(-ycur)*tan(phie)) );
		for(int zcur=zfm; zcur<=zto; zcur++)
		  for(int xcur=xfm; xcur<=xto; xcur++) {
			int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
			int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
			int tmpz = zcur%zn;				  if(tmpz<-zh) tmpz+=zn;				  if(tmpz>=-zh+zn) tmpz-=zn;

			float ss = sma1(xcur)*sma2(ycur)*sma3(zcur);
			float bb = big1(xcur)*big2(ycur)*big3(zcur);

			wpdata(tmpx, tmpy, tmpz) = O(xcur,ycur,zcur) * bb * sqrtf(1.0-ss*ss);

			float thtcur = atan2(zcur/R3, (-ycur)/R2);
			float phicur = atan2(xcur/R1, (-ycur)/R2);

			float glbpou;
			fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou); //CHECK
			float wtht;
			if(thtcur<thtm) {
			  if(h==0)						wtht = 1;
			  else {
				  	 float l,r;
				  	 fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);
				  	 wtht = l;
			  }
			} else {
			  if(h==nd-1)						wtht = 1;
			  else {
				  float l,r;
				  fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);
				  wtht = r;
			  }
			}
			float wphi;
			if(phicur<phim) {
			  if(f==0)						wphi = 1;
			  else {
				  float l,r;
				  fdct3d_window( (phicur-phis)/(phim-phis), l, r);
				  wphi = l;
			  }
			} else {
			  if(f==nd-1)						wphi = 1;
			  else {
				  float l,r;
				  fdct3d_window( (phicur-phim)/(phie-phim), l, r);
				  wphi = r;
			  }
			}
			float pou = glbpou * wtht * wphi;
			wpdata(tmpx, tmpy, tmpz) *= pou;
		  }
	  } //ycur
	  CpxNumTns tpdata(xn,yn,zn);
	  fdct3d_ifftshift(xn,yn,zn,wpdata,tpdata);
	  fftwf_plan p = NULL;

	  p = fftwf_plan_dft_3d(zn, yn, xn, (fftwf_complex*)tpdata.data(), (fftwf_complex*)tpdata.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
	  fftwf_execute(p);		  //cerr<<"wedge s"<<endl;

	  float sqrtprod = sqrt(float(xn*yn*zn));
	  for(int i=0; i<xn; i++)
		  for(int j=0; j<yn; j++)
			  for(int k=0; k<zn; k++)
				  tpdata(i,j,k) /= sqrtprod;

	  fftwf_destroy_plan(p);
	  cellStruct[s].push_back(vector<int> ());
	  cellStruct[s][wcnt].push_back(xn);
	  cellStruct[s][wcnt].push_back(yn);
	  cellStruct[s][wcnt].push_back(zn);

	  csc[wcnt] = tpdata;
	  wcnt++;
	}
  }//end of face
  //face 5. -z,-x,-y
  for(int g=nd-1; g>=0; g--) {
	for(int f=nd-1; f>=0; f--) {
	  
	  float zs = -R3;		  float ze = -R3/4+(W3/2)/4;
	  float xs = -R1 + (2*f-1)*W1/2;		float xe = -R1 + (2*f+3)*W1/2;
	  float ys = -R2 + (2*g-1)*W2/2;		float ye = -R2 + (2*g+3)*W2/2;

	  int xn = int(ceil(xe-xs));		  int yn = int(ceil(ye-ys));		  int zn = int(ceil(ze-zs));
	  float thts, thtm, thte; //y to x
	  if(f==0) {
		thts = atan2(-1.0, 1.0-1.0/nd);
		thtm = atan2(-1.0+1.0/nd, 1.0);
		thte = atan2(-1.0+3.0/nd, 1.0);
	  } else if(f==nd-1) {
		thts = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);
		thtm = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);
		thte = atan2(1.0, 1.0-1.0/nd);
	  } else {
		thts = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);
		thtm = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);
		thte = atan2(-1.0+(2.0*f+3.0)/nd, 1.0);
	  }
	  float phis, phim, phie; //z to x
	  if(g==0) {
		phis = atan2(-1.0, 1.0-1.0/nd);
		phim = atan2(-1.0+1.0/nd, 1.0);
		phie = atan2(-1.0+3.0/nd, 1.0);
	  } else if(g==nd-1) {
		phis = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);
		phim = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);
		phie = atan2(1.0, 1.0-1.0/nd);
	  } else {
		phis = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);
		phim = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);
		phie = atan2(-1.0+(2.0*g+3.0)/nd, 1.0);
	  }
	  int xh = xn/2;		  int yh = yn/2;		  int zh = zn/2;
	  float R13 = R1/R3;		  float R23 = R2/R3;	  //float R13 = float(F1)/float(F3);		  float R23 = float(F2)/float(F3);
	  CpxOffTns wpdata(xn,yn,zn);
	  for(int zcur=(int)ceil(zs); zcur<ze; zcur++) {
		int xfm = (int)ceil( max(-R1, R13*(-zcur)*tan(thts)) );
		int xto = (int)floor( min(R1, R13*(-zcur)*tan(thte)) );
		int yfm = (int)ceil( max(-R2, R23*(-zcur)*tan(phis)) );
		int yto = (int)floor( min(R2, R23*(-zcur)*tan(phie)) );
		for(int xcur=xfm; xcur<=xto; xcur++)
		  for(int ycur=yfm; ycur<=yto; ycur++)
		  {
			int tmpx = xcur%xn;				  if(tmpx<-xh) tmpx+=xn;				  if(tmpx>=-xh+xn) tmpx-=xn;
			int tmpy = ycur%yn;				  if(tmpy<-yh) tmpy+=yn;				  if(tmpy>=-yh+yn) tmpy-=yn;
			int tmpz = zcur%zn;				  if(tmpz<-zh) tmpz+=zn;				  if(tmpz>=-zh+zn) tmpz-=zn;
			float ss = sma1(xcur)*sma2(ycur)*sma3(zcur);
			float bb = big1(xcur)*big2(ycur)*big3(zcur);

			wpdata(tmpx, tmpy, tmpz) = O(xcur,ycur,zcur) * bb * sqrtf(1.0-ss*ss);
			float thtcur = atan2(xcur/R1, (-zcur)/R3);
			float phicur = atan2(ycur/R2, (-zcur)/R3);
			float glbpou;
			fdct3d_globalpou(thtcur, phicur, M_PI/4-atan2(1.0-1.0/nd, 1.0), glbpou);
			float wtht;
			if(thtcur<thtm) {
			  if(f==0)						wtht = 1;
			  else {
				  float l,r;
				  fdct3d_window( (thtcur-thts)/(thtm-thts), l, r);
				  wtht = l;
			  }
			} else {
			  if(f==nd-1)						wtht = 1;
			  else {
				  float l,r;
				  fdct3d_window( (thtcur-thtm)/(thte-thtm), l, r);
				  wtht = r;
			  }
			}
			float wphi;
			if(phicur<phim) {
			  if(g==0)						wphi = 1;
			  else {
				  float l,r;
				  fdct3d_window( (phicur-phis)/(phim-phis), l, r);
				  wphi = l;
			  }
			} else {
			  if(g==nd-1)						wphi = 1;
			  else {
				  float l,r;
				  fdct3d_window( (phicur-phim)/(phie-phim), l, r);
				  wphi = r;
			  }
			}
			float pou = glbpou * wtht * wphi;
			wpdata(tmpx, tmpy, tmpz) *= pou;
		  }
	  }//zcur

	  CpxNumTns tpdata(xn,yn,zn);
	  fdct3d_ifftshift(xn,yn,zn,wpdata,tpdata);
	  fftwf_plan p = NULL;
	  p = fftwf_plan_dft_3d(zn, yn, xn, (fftwf_complex*)tpdata.data(), (fftwf_complex*)tpdata.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
	  fftwf_execute(p);		  //cerr<<"wedge s"<<endl;

	  float sqrtprod = sqrt(float(xn*yn*zn));
	  for(int i=0; i<xn; i++)
		  for(int j=0; j<yn; j++)
			  for(int k=0; k<zn; k++)
				  tpdata(i,j,k) /= sqrtprod;

	  fftwf_destroy_plan(p);

	  cellStruct[s].push_back(vector<int> ());
	  cellStruct[s][wcnt].push_back(xn);
	  cellStruct[s][wcnt].push_back(yn);
	  cellStruct[s][wcnt].push_back(zn);

	  csc[wcnt] = tpdata;
	  wcnt++;
	}
  }//end of face

  assert(wcnt==nd*nd*nf);
  return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------
int fdct3d_forward_wavelet(float L1, float L2, float L3, int s, CpxOffTns& O, vector< vector<CpxNumTns> >& C, std::vector<std::vector<std::vector<int> > > &cellStruct)
{
  vector<CpxNumTns>& csc = C[s];
  cellStruct.push_back(vector<vector<int> > ());
  cellStruct[s].push_back(vector<int> ());
  csc.resize(1);
  
  L1 = L1/2;  L2 = L2/2;  L3 = L3/2;
  int S1, S2, S3;
  int F1, F2, F3;
  float R1, R2, R3;
  fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);

  DblOffVec big1(S1);  fdct3d_lowpass(L1, big1);
  DblOffVec big2(S2);  fdct3d_lowpass(L2, big2);
  DblOffVec big3(S3);  fdct3d_lowpass(L3, big3);
  
  for(int i=-S1/2; i<-S1/2+S1; i++)
	for(int j=-S2/2; j<-S2/2+S2; j++)
	  for(int k=-S3/2; k<-S3/2+S3; k++)
	  {
		float pou = big1(i)*big2(j)*big3(k);
		O(i,j,k) = O(i,j,k) * sqrt(1-pou*pou);
	  }
  //cerr<<energy(O)<<endl;

  int N1 = O.m();  int N2 = O.n();  int N3 = O.p();
  CpxNumTns T(N1,N2,N3);

  fdct3d_ifftshift(N1,N2,N3,O,T);

  fftwf_plan p = fftwf_plan_dft_3d(N3,N2,N1, (fftwf_complex*)T.data(), (fftwf_complex*)T.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwf_execute(p);
  fftwf_destroy_plan(p);

  float sqrtprod = sqrt(float(N1*N2*N3));
  
  for(int i=0; i<N1; i++)
	  for(int j=0; j<N2; j++)
		  for(int k=0; k<N3; k++)
			  T(i,j,k) /= sqrtprod;

  csc[0] = T;

  cellStruct[s][0].push_back(N1);
  cellStruct[s][0].push_back(N2);
  cellStruct[s][0].push_back(N3);

  return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------
int fdct3d_forward_center(float L1, float L2, float L3, int s, CpxOffTns& O, vector< vector<CpxNumTns> >& C, std::vector<std::vector<std::vector<int> > > &cellStruct)
{
  vector<CpxNumTns>& csc = C[s];
  cellStruct.push_back(vector<vector<int> > ());
  cellStruct[0].push_back(vector<int> ());

  csc.resize(1);
  
  int S1, S2, S3;
  int F1, F2, F3;
  float R1, R2, R3;
  fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);

  DblOffVec big1(S1);  fdct3d_lowpass(L1, big1);
  DblOffVec big2(S2);  fdct3d_lowpass(L2, big2);
  DblOffVec big3(S3);  fdct3d_lowpass(L3, big3);

  CpxOffTns A(S1,S2,S3);

  for(int i=-S1/2; i<-S1/2+S1; i++)
	for(int j=-S2/2; j<-S2/2+S2; j++)
	  for(int k=-S3/2; k<-S3/2+S3; k++)
		A(i,j,k) = O(i,j,k) * (big1(i)*big2(j)*big3(k));
  //cerr<<energy(A)<<endl;
  
  CpxNumTns T(S1,S2,S3);
  fdct3d_ifftshift(S1,S2,S3,A,T);

  fftwf_plan p = fftwf_plan_dft_3d(S3,S2,S1, (fftwf_complex*)T.data(), (fftwf_complex*)T.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwf_execute(p);
  fftwf_destroy_plan(p);

  float sqrtprod = sqrtf(float(S1*S2*S3));
  for(int i=0; i<S1; i++)
	  for(int j=0; j<S2; j++)
		  for(int k=0; k<S3; k++)
			  T(i,j,k) /= sqrtprod;

 //// CSC is the 1-D vector of class numtns ***********************************************
  csc[0] = T;
  cellStruct[0][0].push_back(S1);
  cellStruct[0][0].push_back(S2);
  cellStruct[0][0].push_back(S3);
  
  return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------
int fdct3d_forward(int N1, int N2, int N3, int nbscales, int nbdstz_coarse, int ac, CpxNumTns& X, vector< vector<CpxNumTns> >& C, std::vector<std::vector<std::vector<int> > > &cellStruct)
{
  /// Make a duplicate of class X to Class T
  CpxNumTns T(X);
  fftwf_plan p = fftwf_plan_dft_3d(N3, N2, N1, (fftwf_complex*)T.data(), (fftwf_complex*)T.data(), FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(p);

  float sqrtprod = sqrt(float(N1*N2*N3));

  for(int i=0; i<N1; i++)
	  for(int j=0; j<N2; j++)
		  for(int k=0; k<N3; k++)
		  {
		      /// T(i,j,k) ==> T.data(i + j*N1 + k*N2*N1) ===> Converting to 3-D to 1-D data
			  T(i,j,k) /= sqrtprod;
		  }

  fftwf_destroy_plan(p);

  /// Generate a class with 6 variables and a long one dimensional pointer with the size of (N1* N2* N3) ***************
  CpxOffTns F(N1, N2, N3);
  fdct3d_fftshift(N1, N2, N3, T, F);
  
  //expand if necessary
  CpxOffTns O;

  if(ac==1)
  {
		float L1 = 4.0*N1/3.0;
		float L2 = 4.0*N2/3.0;
		float L3 = 4.0*N3/3.0;

		int S1, S2, S3;
		int F1, F2, F3;
		float R1, R2, R3;
		fdct3d_rangecompute(L1, L2, L3, S1, S2, S3, F1, F2, F3, R1, R2, R3);

		IntOffVec t1(S1);
		for(int i=-F1; i<-F1+S1; i++)
			if(i<-N1/2)
				t1(i) = i+int(N1);
			else if(i>(N1-1)/2)
				t1(i) = i-int(N1);
			else t1(i) = i;

		IntOffVec t2(S2);
		for(int i=-F2; i<-F2+S2; i++)
			if(i<-N2/2)
				t2(i) = i+int(N2);
			else if(i>(N2-1)/2)
				t2(i) = i-int(N2);
			else t2(i) = i;

		IntOffVec t3(S3);
		for(int i=-F3; i<-F3+S3; i++)
			if(i<-N3/2)
				t3(i) = i+int(N3);
			else if(i>(N3-1)/2)
				t3(i) = i-int(N3);
			else t3(i) = i;

		O.resize(S1, S2, S3);
		for(int i=-F1; i<-F1+S1; i++)
			for(int j=-F2; j<-F2+S2; j++)
				for(int k=-F3; k<-F3+S3; k++)
					O(i,j,k) = F(t1(i), t2(j), t3(k));
  }
  else
  {
	   O = F;
  }
  
  int L = nbscales;
  /// Resize the outer vector of C variable  **********************************************8
  C.resize(L);

  if(ac==1) {
	{
	  int s = 0;
	  float L1 = 4.0*N1/3.0 / pow2(L-1-s);	 float L2 = 4.0*N2/3.0 / pow2(L-1-s);	 float L3 = 4.0*N3/3.0 / pow2(L-1-s);
	  fdct3d_forward_center(L1,L2,L3,s, O, C, cellStruct);
	}

	for(int s=1; s<L; s++) {
	  float L1 = 4.0*N1/3.0 / pow2(L-1-s);	 float L2 = 4.0*N2/3.0 / pow2(L-1-s);	 float L3 = 4.0*N3/3.0 / pow2(L-1-s);
	  int nd = nbdstz_coarse * pow2(s/2);
	  fdct3d_forward_angles(L1,L2,L3,s, nd, O, C, cellStruct);
	}

  } else {
	{
	  int s = 0;
	  float L1 = 4.0*N1/3.0 / pow2(L-1-s);	 float L2 = 4.0*N2/3.0 / pow2(L-1-s);	 float L3 = 4.0*N3/3.0 / pow2(L-1-s);
	  fdct3d_forward_center(L1,L2,L3,s, O, C, cellStruct);
	}
	for(int s=1; s<L-1; s++) {
	  float L1 = 4.0*N1/3.0 / pow2(L-1-s);
	  float L2 = 4.0*N2/3.0 / pow2(L-1-s);
	  float L3 = 4.0*N3/3.0 / pow2(L-1-s);
	  int nd = nbdstz_coarse * pow2(s/2);
	  fdct3d_forward_angles(L1,L2,L3,s, nd, O, C, cellStruct);
	}
	{
	  int s = L-1;
	  float L1 = 4.0*N1/3.0 / pow2(L-1-s);
	  float L2 = 4.0*N2/3.0 / pow2(L-1-s);
	  float L3 = 4.0*N3/3.0 / pow2(L-1-s);
	  fdct3d_forward_wavelet(L1,L2,L3,s, O, C, cellStruct);
	}
  }

  O.~OffTns();
  F.~OffTns();

  return 0;
}
