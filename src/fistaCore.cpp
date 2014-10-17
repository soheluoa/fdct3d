/*
 *  =================================================================================
 *  Name        : fistaCore.cpp
 *  Author      : Sohel Bhuiyan
 *  Version     : 1.0
 *  Purpose     : Core cpp file to apply FISTA.
 *  Date        : September 12, 2014
 *  Affiliation : University of Alberta, Physics department, (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */

#include <math.h>
#include <stdlib.h>

#include "fistaParams.hpp"
#include "fistaCore.hpp"
#include "PARAMS.hpp"
#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

#define PI 3.1415926535897932

/******** Core method to execute FISTA algorithm ****************************************/
void fistaCore::Reconstruct(fistaParams* fsp, PARAMS* params, fistaCore* fsc)
{
    int chkDataFile=0;

	/**** Read a binary file (Seismic data) ********************************************/
    chkDataFile = readBinFile(fsp->inDataFileName, fsp->n1, fsp->n2, fsp->n3, params->inData);

    // If the binary file is read successfully *****************************************/
	if(chkDataFile == 1)
	{
		double tmpt, temp, normType, t=1;
		CpxNumTns obsData(fsp->n3,fsp->n2,fsp->n1);					// y
		CpxNumTns tmpData(fsp->n3,fsp->n2,fsp->n1);                 // Hx
		CpxNumTns diffData(fsp->n3,fsp->n2,fsp->n1);                // y-Hx
		CpxNumTns reconData(fsp->n3,fsp->n2,fsp->n1);               // recon data
		std::vector< std::vector<CpxNumTns> > curvCoeff;            // x
		std::vector< std::vector<CpxNumTns> > tempCurvCoeff;        // Yk
		std::vector< std::vector<CpxNumTns> > addCurvCoeff;         // temp_Yk_mat
		std::vector< std::vector<CpxNumTns> > iterTempCurvCoeff;    // tempx
		std::vector< std::vector<CpxNumTns> > iterDiffTempCurvCoeff;// diff coeff

		computeObsData(fsp->sampleMatFileName, fsp->n1, fsp->n2, fsp->n3, params->inData, params->samplMat, params->obsData);

        // Save seismic data as complex number *****************************************
		for(int i=0; i<fsp->n3; i++)
		   for(int j=0; j<fsp->n2; j++)
			   for(int k=0; k<fsp->n1; k++)
				   obsData(i,j,k) =  cpx(params->obsData[i][j][k],0);

		// Initialise the parameters for 3D curvelet transformation ********************
		fdct3d_param(fsp->n3,fsp->n2,fsp->n1,params->nbscales,params->nbdstz_coarse,params->ac,params->fxs,params->fys,params->fzs, params->nxs,params->nys,params->nzs);

		// Apply forward curvelet operator *********************************************
		fdct3d_forward(fsp->n3,fsp->n2,fsp->n1,params->nbscales,params->nbdstz_coarse, params->ac, obsData, curvCoeff, params->cellStruct);

		// Traversing the curvelet coefficients and assign zero ************************
		for (unsigned int s1=0; s1<params->cellStruct.size(); s1++)
		{
			tempCurvCoeff.push_back(vector<CpxNumTns> ());
			for (unsigned int s2=0; s2<params->cellStruct[s1].size(); s2++)
			{
				CpxNumTns tmpCoeff(params->cellStruct[s1][s2][0],params->cellStruct[s1][s2][1],params->cellStruct[s1][s2][2]);
				tmpCoeff = curvCoeff[s1][s2];
				for (int i=0; i<params->cellStruct[s1][s2][0]; i++)
					for (int j=0; j<params->cellStruct[s1][s2][1]; j++)
						for (int k=0; k<params->cellStruct[s1][s2][2]; k++)
							tmpCoeff(i,j,k) =cpx(0.0, 0.0);
				curvCoeff[s1][s2] = tmpCoeff;
				tempCurvCoeff[s1].push_back(tmpCoeff);
				tmpCoeff.~NumTns();
			}
		}

		// Starting iterations to converge the FISTA algorithm ************************
		for(int itrNum=0; itrNum<fsp->maxItr; itrNum++)
		{
			// Traversing the curvelet coefficients and Yk = x ************************
			for (unsigned int s1=0; s1<params->cellStruct.size(); s1++)
			{
				iterTempCurvCoeff.push_back(vector<CpxNumTns> ());
				for (unsigned int s2=0; s2<params->cellStruct[s1].size(); s2++)
				{
					CpxNumTns tmpCoeff(params->cellStruct[s1][s2][0],params->cellStruct[s1][s2][1],params->cellStruct[s1][s2][2]);
					tmpCoeff = curvCoeff[s1][s2];
					iterTempCurvCoeff[s1].push_back(tmpCoeff);
					tmpCoeff.~NumTns();
				}
			}

			// Apply adjoint or transposed curvelet operator; y = Ax; input: tempCurvCoeff, output: tmpData
			fdct3d_inverse(fsp->n3,fsp->n2,fsp->n1,params->nbscales,params->nbdstz_coarse, params->ac, tempCurvCoeff, tmpData);
			dotProduct(fsp->n1, fsp->n2, fsp->n3, params->samplMat, tmpData);

			for(int i=0; i<fsp->n3; i++)
				for(int j=0; j<fsp->n2; j++)
					for(int k=0; k<fsp->n1; k++)
						diffData(i,j,k) = (real(obsData(i,j,k)) - real(tmpData(i,j,k)));

			params->cellStruct.clear();

			// x = A'y; y = y - Hx; input: params->diffData, output: params->tempDiffCurvCoeff
			fdct3d_forward(fsp->n3,fsp->n2,fsp->n1,params->nbscales,params->nbdstz_coarse, params->ac, diffData, iterDiffTempCurvCoeff, params->cellStruct);

			// Traversing the curvelet coefficients and assign zero ************************
			for (unsigned int s1=0; s1<params->cellStruct.size(); s1++)
			{
				addCurvCoeff.push_back(vector<CpxNumTns> ());
				for (unsigned int s2=0; s2<params->cellStruct[s1].size(); s2++)
				{
					CpxNumTns tmpCoeff1(params->cellStruct[s1][s2][0],params->cellStruct[s1][s2][1],params->cellStruct[s1][s2][2]);
					CpxNumTns tmpCoeff2(params->cellStruct[s1][s2][0],params->cellStruct[s1][s2][1],params->cellStruct[s1][s2][2]);
					CpxNumTns tmpCoeff(params->cellStruct[s1][s2][0],params->cellStruct[s1][s2][1],params->cellStruct[s1][s2][2]);
					tmpCoeff1 = iterDiffTempCurvCoeff[s1][s2];
					tmpCoeff2 = tempCurvCoeff[s1][s2];
					for (int i=0; i<params->cellStruct[s1][s2][0]; i++)
						for (int j=0; j<params->cellStruct[s1][s2][1]; j++)
							for (int k=0; k<params->cellStruct[s1][s2][2]; k++)
								tmpCoeff(i,j,k) = cpx(real(tmpCoeff1(i,j,k)) + (real(tmpCoeff2(i,j,k))/fsp->alpha), 0.0);
					addCurvCoeff[s1].push_back(tmpCoeff);
					tmpCoeff1.~NumTns();
					tmpCoeff2.~NumTns();
					tmpCoeff.~NumTns();
				}
			}

			wthresh(fsp->thresh, addCurvCoeff, curvCoeff, params->cellStruct);
			addCurvCoeff.clear();

			tmpt = t;
			t = (1+sqrt(1+4*powf(t,2)))/2;
			temp = ((tmpt-1)/t);

			// Traversing the curvelet coefficients and assign zero ************************
			for (unsigned int s1=0; s1<params->cellStruct.size(); s1++)
			{
				for (unsigned int s2=0; s2<params->cellStruct[s1].size(); s2++)
				{
					CpxNumTns tmpCoeff1(params->cellStruct[s1][s2][0],params->cellStruct[s1][s2][1],params->cellStruct[s1][s2][2]);
					CpxNumTns tmpCoeff2(params->cellStruct[s1][s2][0],params->cellStruct[s1][s2][1],params->cellStruct[s1][s2][2]);
					CpxNumTns tmpCoeff(params->cellStruct[s1][s2][0],params->cellStruct[s1][s2][1],params->cellStruct[s1][s2][2]);

					tmpCoeff1 = curvCoeff[s1][s2];
					tmpCoeff2 = iterTempCurvCoeff[s1][s2];
					for (int i=0; i<params->cellStruct[s1][s2][0]; i++)
						for (int j=0; j<params->cellStruct[s1][s2][1]; j++)
							for (int k=0; k<params->cellStruct[s1][s2][2]; k++)
								tmpCoeff(i,j,k) =  cpx(real(tmpCoeff1(i,j,k)) + temp * (real(tmpCoeff1(i,j,k)) - real(tmpCoeff2(i,j,k))), 0.0);
					tempCurvCoeff[s1][s2] = tmpCoeff;
					tmpCoeff1.~NumTns();
					tmpCoeff2.~NumTns();
					tmpCoeff.~NumTns();
				}
			}
			iterTempCurvCoeff.clear();
			std::cout << "Iteration completed" << " "<< itrNum <<std::endl;
		}

		// y = Ax; input: params->tempCurvCoeff, output: params->reconData
		fdct3d_inverse(fsp->n3,fsp->n2,fsp->n1,params->nbscales,params->nbdstz_coarse, params->ac, tempCurvCoeff, reconData);

		// Writing the binary file for the reconstructed data
		writeBinFile(fsp->outDataFileName, fsp->n1, fsp->n2, fsp->n3, reconData);

		obsData.~NumTns();					// y
		tmpData.~NumTns();                   // Hx
		diffData.~NumTns();                 // y-Hx
		reconData.~NumTns();                // recon data
		curvCoeff.clear();            		// x
		tempCurvCoeff.clear();        		// Yk
		addCurvCoeff.clear();         		// temp_Yk_mat
		iterTempCurvCoeff.clear();    		// tempx
		iterDiffTempCurvCoeff.clear();		// diff coeff

		system("cd /home/entropy/workspace/fdct3d/src/");
		system("./reconCube.sh");
	}
	else
		std::cout << "Please check whether the path of binary files are correct or not. " <<std::endl;
}

/*  3D complex float dynamic memory allocation with the size of ( n1 X n2 x n3 ) ***************/
inline int fistaCore::readBinFile(std::string fileName, int n1, int n2, int n3, vector<vector<vector<double > > > &data)
{

	int fsize, chkFile = 0;
	FILE *fp;
	fp = fopen(fileName.c_str(), "rb");
	double *temp, maxVal;

    if(!fp)
    	std::cout << "File could not opened. " <<std::endl;
    else
    {
    	chkFile = 1;
    	fsize = n1 * n2 * n3;
    	temp = (double *) malloc(sizeof(double) * fsize);
    	fread(temp,sizeof(double), fsize, fp);
    	double max = *std::max_element(temp,temp+fsize);

        // Very important: vector data stores the 3D seismic traces as follows:
    	// k-->Time samples, j-->(n2) in pscube, i-->n3 in pscube
		for(int i=0; i<n3; i++)
		{
			data.push_back(vector<vector<double> > ());
			for(int j=0; j<n2; j++)
			{
			   data[i].push_back(vector<double> ());
			   for(int k=0; k<n1; k++)
				   data[i][j].push_back(temp[i*n2*n1+j*n1+k]*max);
			}
		 }
    }
	return chkFile;
}


inline int fistaCore::dotProduct(int n1, int n2, int n3, std::vector<std::vector<int> > &sampleMat, CpxNumTns &data )
{
	for(int i=0; i<n3; i++)
	{
	   for(int j=0; j<n2; j++)
	   {
		   if(sampleMat[i][j] == 0)
		   {
			   for( int k=0; k<n1; k++)
				   data(i,j,k) = 0.0;
		   }
		   else
		   {
			   for( int k=0; k<n1; k++)
			   	   data(i,j,k) = cpx(real(data(i,j,k)), 0.0);
		   }
	   }
	}

	return 0;
}


/*  3D complex float observed data (Orig_data .* sampling ) with the size of ( n1 X n2 x n3 ) ***************/
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

/* Compute different types of norms  *****************************************************/
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

/******* Thresholding operator for the curvelet coefficients ************************************************/
inline void fistaCore::wthresh(double thresh, std::vector<std::vector<CpxNumTns > >&x, std::vector<std::vector<CpxNumTns > >&y, std::vector<std::vector<std::vector<int> > > &cellStruct )
{
    double temp;
    int sgn_data;

    // Traversing the curvelet coefficients and assign zero ************************
	for (unsigned int s1=0; s1<cellStruct.size(); s1++)
	{
		for (unsigned int s2=0; s2<cellStruct[s1].size(); s2++)
		{
			CpxNumTns tmpCoeff(cellStruct[s1][s2][0],cellStruct[s1][s2][1],cellStruct[s1][s2][2]);
			tmpCoeff = x[s1][s2];
			for (int i=0; i<cellStruct[s1][s2][0]; i++)
				for (int j=0; j<cellStruct[s1][s2][1]; j++)
					for (int k=0; k<cellStruct[s1][s2][2]; k++)
					{
						 temp = (fabs(real(tmpCoeff(i,j,k)))-thresh);
						 temp = (temp+abs(temp))/2;
						 sgn_data = ((real(tmpCoeff(i,j,k)) > 0) ? 1 : (real(tmpCoeff(i,j,k)) < 0) ? -1 : 0);
						 tmpCoeff(i,j,k) = cpx(sgn_data * temp, 0.0);
					}
			y[s1][s2] = tmpCoeff;
			tmpCoeff.~NumTns();
		}
	}
}

/*  Write the binary file of the reconstructed seismic data with the size of ( n1 X n2 X n3 )***************/
inline void fistaCore::writeBinFile(std::string fileName,int n1, int n2, int n3, CpxNumTns &data)
{
	int fsize;
	FILE *fp;
	fp = fopen(fileName.c_str(), "rb");
	double *temp;

	fsize = n1 * n2* n3;
	temp = (double *) malloc(sizeof(double) * fsize);

	for(int i=0; i<n3; i++)
		for(int j=0; j<n2; j++)
			for(int k=0; k<n1; k++)
				temp[i*n2*n1+j*n1+k] = real(data(i,j,k));

	fp = fopen(fileName.c_str(), "w");
	fwrite(temp,sizeof(float),fsize,fp);
	fclose(fp);
}
