/*
 * fistaCore.cpp
 *
 *  Created on: Sep 14, 2014
 *      Author: entropy
 */

#include <cmath>
#include <math.h>
#include <algorithm>
using::std::sort;
#include <stdlib.h>
#include <vector>
using std::vector;

#include <iostream>
using std::cout;
using std::ofstream;
using std::endl;

#include <fstream>
#include <string>

#include <complex.h>
#include <fftw3.h>

#include "fistaParams.hpp"
#include "fistaCore.hpp"
#include "PARAMS.hpp"

#define PI 3.1415926535897932

/* Core method to generate seismic image using PSPI algorithm *****************/
void fistaCore::Execute(fistaParams* fsp, PARAMS* params, fistaCore* fsc)
{
    int chkDataFile=0, normType, t=1;

	/* Read a binary file (Seismic data) **************************************/
    chkDataFile = readBinFile(fsp->inDataFileName, fsp->n1, fsp->n2, fsp->n3, params->inData);
    normType = 1;

	if(chkDataFile == 1)
	{
		double tmpt, temp;

		computeObsData(fsp->sampleMatFileName, fsp->n1, fsp->n2, fsp->n3, params->inData, params->samplMat, params->obsData);

		/*Write the reconstructed seismic data ********************************/
		writeBinFile(fsp->outDataFileName, fsp->n1, fsp->n2, fsp->n3, params->inData);

		writeBinFile(fsp->outDataFileName, fsp->n1, fsp->n2, fsp->n3, params->obsData);

		/* x = A'y; input: params->obsData, output: params->curvCoeff */
		aprod(fsp, params, -1);

		for(int i=0; i<params->Nw; i++)
			for(int j=0; j<params->Kx; j++)
				for(int k=0; k<params->Ky; k++)
				{
					params->curvCoeff[i][j][k] = 0*params->curvCoeff[i][j][k];
					params->tempCurvCoeff[i][j][k] = params->curvCoeff[i][j][k];
				}

		for(int i=0; i<fsp->maxItr; i++)
		{

			for(int i=0; i<params->Nw; i++)
				for(int j=0; j<params->Kx; j++)
					for(int k=0; k<params->Ky; k++)
						params->tempModel[i][j][k] = params->curvCoeff[i][j][k];

			/* y = Ax; Hx; input: params->tempCurvCoeff, output: params->tempData */
			aprod(fsp, params, 1);

			for(int i=0; i<fsp->n3; i++)
				for(int j=0; j<fsp->n2; j++)
					for(int k=0; k<fsp->n1; k++)
						params->tempData[i][j][k] = (params->obsData[i][j][k] - params->tempData[i][j][k]);

			/* x = A'y; y = y - Hx; input: params->tempData, output: params->tempDiffCurvCoeff */
			aprod(fsp, params, -1);

			for(int i=0; i<params->Nw; i++)
				for(int j=0; j<params->Kx; j++)
					for(int k=0; k<params->Ky; k++)
						params->newModel[i][j][k] = params->tempCurvCoeff[i][j][k] + params->tempDiffCurvCoeff[i][j][k]/fsp->alpha;

			wthresh(fsp->n1, fsp->n2, fsp->n3,  params->newModel, fsp->thresh, params->curvCoeff);

			tmpt = t;

			t = (1+sqrt(1+4*powf(t,2)))/2;

			temp = ((tmpt-1)/t);

			for(int i=0; i<params->Nw; i++)
				for(int j=0; j<params->Kx; j++)
					for(int k=0; k<params->Ky; k++)
						params->tempCurvCoeff[i][j][k] = params->curvCoeff[i][j][k] + temp* (params->curvCoeff[i][j][k] - params->tempModel[i][j][k]);

		}

		/* y = Ax; input: params->tempCurvCoeff, output: params->reconData */
		aprod(fsp, params, 1);

		/* Writing the binary file for the reconstructed data */
		writeBinFile(fsp->outDataFileName, fsp->n1, fsp->n2, fsp->n3, params->reconData);
		system("reconCube.sh");

		//normVal = norm(params->samplMat, normType);
	}
	else
		std::cout << "Please check whether the path of binary files are correct or not. " <<std::endl;
}

/*  2D complex float dynamic memory allocation with n1 X n2 size ***************/
inline int fistaCore::readBinFile(std::string fileName, int n1, int n2, int n3, vector<vector<vector<double > > > &data)
{

	int fsize, chkFile = 0;
	FILE *fp;
	std::cout << fileName << std::endl;
	fp = fopen(fileName.c_str(), "rb");
	double *temp;

    if(!fp)
    	std::cout << "File could not opened. " <<std::endl;
    else
    {
    	chkFile = 1;
    	fsize = n1 * n2 * n3;
    	temp = (double *) malloc(sizeof(double) * fsize);
    	fread(temp,sizeof(double), fsize, fp);

        // Very important: vector data stores the 3D seismic traces as follows:
    	// k-->Time samples, j-->(n2) in pscube, i-->n3 in pscube
		for(int i=0; i<n3; i++)
		{
			data.push_back(vector<vector<double> > ());
			for(int j=0; j<n2; j++)
			{
			   data[i].push_back(vector<double> ());
			   for(int k=0; k<n1; k++)
			   {
				   data[i][j].push_back(temp[i*n2*n1+j*n1+k]);
			   }
			}
		 }
    }
	return chkFile;
}



inline int fistaCore::computeObsData(std::string fileName, int n1, int n2, int n3, std::vector<std::vector<std::vector<double> > > &data, std::vector<std::vector<int> > &sampleMat, std::vector<std::vector<std::vector<double> > > &obsData )
{
	int chkDone = 0, fsize, *temp;
	FILE *fp;
	fp = fopen(fileName.c_str(), "rb");

    if(!fp)
    	std::cout << "File could not opened. " <<std::endl;
    else
    {
    	chkDone = 1;
    	fsize = n2 * n3;
    	temp = (int *) malloc(sizeof(int) * fsize);
    	fread(temp,sizeof(int), fsize, fp);

    	for(int i=0; i<n3; i++)
		{
			obsData.push_back(vector<vector<double> > ());
			for(int j=0; j<n2; j++)
			{
			   obsData[i].push_back(vector<double> ());
			   for(int k=0; k<n1; k++)
				   obsData[i][j].push_back(data[i][j][k]);
			}
		}

    	for(int i=0; i<n3; i++)
		{
           sampleMat.push_back(vector<int> ());
		   for(int j=0; j<n2; j++)
		   {
			   sampleMat[i].push_back(temp[i*n2+j]);
		   	   if(sampleMat[i][j] == 0)
		   	   {
		   		   for( int k=0; k<n1; k++)
		   			   obsData[i][j][k] = 0.0;
		   	   }
		   }
		}
    }
	return chkDone;
}


inline double fistaCore::norm(vector<vector<int> > &x, int normType)
{
	double normVal = 0.0;
	int fsize, *temp;

	if(normType == 2)
	{
		for(unsigned int i=0; i<x.size(); i++)
		{
			for(unsigned int j=0; j<x[0].size(); j++)
			{
			   //for(unsigned int k=0; k<x[0][0].size(); k++)
				   normVal += x[i][j]*x[i][j];
			}
		}
		normVal = sqrt (normVal);
	}

	if(normType == 1)
	{
		for(unsigned int i=0; i<x.size(); i++)
		{
			for(unsigned int j=0; j<x[0].size(); j++)
			{
			   //for(unsigned int k=0; k<x[0][0].size(); k++)
				   normVal += abs(x[i][j]);
			}
		}
	}

	// 'Inf' norm type /////
	if(normType == 3)
	{
		fsize = x.size() * x[0].size();
		temp = (int *) malloc(sizeof(int) * fsize);

		for(unsigned int i=0; i<x.size(); i++)
		{
			for(unsigned int j=0; j<x[0].size(); j++)
			{
			   //for(unsigned int k=0; k<x[0][0].size(); k++)
				   temp[i*x[0].size()+j] = (x[i][j]);
			}
		}
		normVal = *std::max_element(temp, temp+fsize);
		normVal = sqrt (normVal);
	}
	return normVal;
}


/* Very hard task to complete; implementing 3D curvelet transformation  */

inline void fistaCore::aprod(fistaParams* fsp, PARAMS* params, int operatorType)
{
	/* Forward operator, Data = A * model, A = sampling * inverseCurvelet */
	if(operatorType == 1)
		int p =0;
		//ifdct_wrapping_sym(&fsp, &params);

	/* Reverse operator, Model = A' * data, A' = Curvelet * sampling'  */
	else if(operatorType == -1)
		int p = 1;
		//fdct_wrapping_sym(fsp->obsData , &params);
}


//inline void spglCore::project(spglParams* fsp, PARAMS* params, int )
//{
//	/* Forward operator, Data = A * model, A = sampling * inverseCurvelet */
//	if(operatorType == 1)
//		int p =0;
//		//ifdct_wrapping_sym(&fsp, &params);
//
//	/* Reverse operator, Model = A' * data, A' = Curvelet * sampling'  */
//	else if(operatorType == -1)
//		int p = 1;
//		//fdct_wrapping_sym(fsp->obsData , &params);
//}


//void spglCore::spgLineCurvy(double f, std::vector<std::vector<std::vector<double > > >x, std::vector<std::vector<std::vector<double > > >dx, double gtd, double fMax, std::vector<std::vector<std::vector<double > > >data)
//{
//	double gamma=1E-04, gts;
//	int maxIts=6, step=1, sNorm=0, scale=1, nSafe=0, iter = 0, err;
//
//    while(1)
//    {
//    	project(spglParams* fsp, PARAMS* params, int );
//    	newMisfit = aprod(&fsp, &params, 1);
//
//    	fNew = misfit
//    	gts = scale * ;
//
//    	if(gts >=0)
//    	{
//    		err = 1;
//    		break;
//    	}
//
//    	if(fNew < fMax + gamma *step*gts)
//    	{
//    		err = 0;
//    		break;
//    	}
//    	else if(iter > maxIts)
//    	{
//    		err = 1;
//    		break;
//    	}
//    	iter++;
//    	step/=2;
//
//    	//    //
//    	sNormOld = sNorm;
//    	if(abs(sNorm - sNormOld) < 1E-06*sNorm)
//    	{
//    		gNorm = norm(g)/sqrt(spg->n1 * spg->n2, spg->n3);
//    		scale = sNorm/gNorm/(pow(2,nSafe));
//    		nSafe++;
//    	}
//    }
//}


inline void fistaCore::wthresh(int n1, int n2, int n3, vector<vector<vector<double> > > &data, double thresh, std::vector<std::vector<std::vector<double > > >&x )
{
    double temp;
    int sgn_data;

	for(int i=0; i<n3; i++)
			for(int j=0; j<n2; j++)
				for(int k=0; k<n1; k++)
				{
					temp = fabs(data[i][j][k]) - thresh;
					temp = (temp+fabs(temp))/2;
					sgn_data = (data[i][j][k] > 0) ? 1 : ((data[i][j][k] < 0) ? -1 : 0);
					x[i][j][k]   = sgn_data*temp;
				}

}



/*  2D float dynamic memory allocation with n1 X n2 size *******************************************/
inline void fistaCore::writeBinFile(std::string fileName,int n1, int n2, int n3, vector<vector<vector<double> > > &data)
{
	int fsize;
	FILE *fp;
	fp = fopen(fileName.c_str(), "rb");
	double *temp;

	std::cout <<fileName <<std::endl;

	fsize = n1 * n2* n3;
	temp = (double *) malloc(sizeof(double) * fsize);

	for(int i=0; i<n3; i++)
		for(int j=0; j<n2; j++)
			for(int k=0; k<n1; k++)
				temp[i*n2*n1+j*n1+k] = data[i][j][k];

	fp = fopen(fileName.c_str(), "w");
	fwrite(temp,sizeof(float),fsize,fp);
	fclose(fp);
}
