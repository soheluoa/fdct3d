/*
 * fistaParams.hpp
 *
 *  Created on: Sep 12, 2014
 *      Author: entropy
 */

#ifndef FISTAPARAMS_HPP_
#define FISTAPARAMS_HPP_

#include<vector>
using std::vector;
#include <iostream>

class fistaParams{

public:
	// These parameters depend on the seismic data and velocity model
	double bpTol;				    			/// Tolerance for basis pursuit solution
    int maxItr;									/// Maximum iterations number
	int nLineTot;	                  		    /// Total no. of linesearch steps.
	double thresh; 								/// Threshold value
	double lambda; 								/// Convergence constant
	double alpha; 								/// Computed using power eigen method
	int n1, n2, n3;
    bool stat;

    std::string inDataFileName;   					/// Absolute path of the In seismic data file
    std::string outDataFileName;   					/// Absolute path of the observed seismic data file
    std::string sampleMatFileName;   				/// Absolute path of the sample data file

	/******** clear configuration variables **********/
	void reset()
	{

	}
};

#endif /* SPGLPARAMS_HPP_ */
