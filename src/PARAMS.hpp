/*
 * PARAMS.hpp
 *
 *  Created on: Sep 15, 2014
 *      Author: entropy
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
	      std::vector<std::vector<std::vector<double > > > reconData;
	      std::vector<std::vector<std::vector<double > > > tempData;
	      std::vector<std::vector<std::vector<double > > > newModel;
	      std::vector< std::vector<double> > fxs,fys,fzs;
	      std::vector< std::vector<int> > nxs,nys,nzs;
	      std::vector<int > cellStruct;

	      // Vector of class NumTns with complex double type ***************************************************


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
	      		reconData.clear();
	      		fxs.clear();
	      		fys.clear();
	      		fzs.clear();
	      		nxs.clear();
	      		nys.clear();
	      		nzs.clear();
	      }
};

#endif /* PARAMS_HPP_ */
