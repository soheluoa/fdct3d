/*
 * PARAMS.hpp
 *
 *  Created on: Sep 15, 2014
 *      Author: entropy
 */

#ifndef PARAMS_HPP_
#define PARAMS_HPP_

#include <vector>

class PARAMS{

public:
	      int Nw;
	      int Kx;
	      int Ky;
	      std::vector<std::vector<int > > samplMat;
	      std::vector<std::vector<int > > cellSize;
	      std::vector<std::vector<std::vector<double > > > curvCoeff;
	      std::vector<std::vector<std::vector<double > > > tempCurvCoeff;
	      std::vector<std::vector<std::vector<double > > > tempModel;
	      std::vector<std::vector<std::vector<double > > > tempDiffCurvCoeff;
	      std::vector<std::vector<std::vector<double > > > inData;
	      std::vector<std::vector<std::vector<double > > > obsData;
	      std::vector<std::vector<std::vector<double > > > reconData;
	      std::vector<std::vector<std::vector<double > > > tempData;
	      std::vector<std::vector<std::vector<double > > > newModel;
	      std::vector<int > cellStruct;

	      void reset()
	      {
	      		samplMat.clear();
	      		cellSize.clear();
	      		curvCoeff.clear();
	      		inData.clear();
	      		obsData.clear();
	      		tempCurvCoeff.clear();
	      		tempModel.clear();
	      		newModel.clear();
	      		tempData.clear();
	      		tempDiffCurvCoeff.clear();
	      		cellStruct.clear();
	      		reconData.clear();
	      }
};

#endif /* PARAMS_HPP_ */
