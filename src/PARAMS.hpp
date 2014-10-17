/*
 *  =================================================================================
 *  Name        : PARAMS.hpp
 *  Author      : Sohel Bhuiyan
 *  Version     : 1.0
 *  Purpose     : Declare the data structure required for the FISTA and 3D curvelet transformation.
 *  Date        : September 12, 2014
 *  Affiliation : University of Alberta, Physics department, (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */

#ifndef PARAMS_HPP_
#define PARAMS_HPP_

#include <vector>
#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

class PARAMS{

public:
	      int Nw;
	      int Kx;
	      int Ky;
	      int nbscales, nbdstz_coarse, ac;
	      std::vector<std::vector<int > > samplMat;
	      std::vector<std::vector<int > > cellSize;
	      std::vector<std::vector<std::vector<double > > > tempCurvCoeff;
	      std::vector<std::vector<std::vector<double > > > tempModel;
	      std::vector<std::vector<std::vector<double > > > tempDiffCurvCoeff;
	      std::vector<std::vector<std::vector<double > > > inData;
	      std::vector<std::vector<std::vector<double > > > obsData;

	      std::vector<std::vector<std::vector<double > > > tempData;
	      std::vector<std::vector<std::vector<double > > > newModel;
	      std::vector< std::vector<double> > fxs,fys,fzs;
	      std::vector< std::vector<int> > nxs,nys,nzs;
	      std::vector<std::vector<std::vector <int > > >cellStruct;

	      void reset()
	      {
	      		samplMat.clear();
	      		cellSize.clear();
	      		inData.clear();
	      		obsData.clear();
	      		tempCurvCoeff.clear();
	      		tempModel.clear();
	      		newModel.clear();
	      		tempData.clear();
	      		tempDiffCurvCoeff.clear();
	      		cellStruct.clear();
	      		fxs.clear();
	      		fys.clear();
	      		fzs.clear();
	      		nxs.clear();
	      		nys.clear();
	      		nzs.clear();
	      }
};

#endif /* PARAMS_HPP_ */
