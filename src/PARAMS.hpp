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
	      int nbscales, nbdstz_coarse, ac;
	      float qFactor;
	      std::vector<std::vector<int > > samplMat;
	      std::vector<double > misFit;
	      std::vector<std::vector<std::vector<float > > > inData;
	      std::vector<std::vector<std::vector<float > > > obsData;
	      std::vector<std::vector<std::vector<float > > > tempData;
	      std::vector< std::vector<float> > fxs,fys,fzs;
	      std::vector< std::vector<int> > nxs,nys,nzs;
	      std::vector<std::vector<std::vector <int > > >cellStruct;

	      void reset()
	      {
	      		samplMat.clear();
	      		inData.clear();
	      		obsData.clear();
	      		tempData.clear();
	      		cellStruct.clear();
	      		fxs.clear();
	      		fys.clear();
	      		fzs.clear();
	      		nxs.clear();
	      		nys.clear();
	      		nzs.clear();
	      		misFit.clear();
	      }
};

#endif /* PARAMS_HPP_ */
