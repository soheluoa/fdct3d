/*
 * fistaCore.cpp
 *
 *  Created on: Sep 12, 2014
 *      Author: entropy
 */
/*  =================================================================================
 *  Name        : FistaCore.hpp
 *  Author      : Sohel Bhuiyan
 *  Version     : 1.0
 *  Purpose     : Base code of this software package. Contains all the basic methods entailed to reconstruct seismic data.
 *  Date        : Feb 18, 2014
 *  Affiliation : University of Alberta, Physics department (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */

#ifndef FISTACORE_HPP_
#define FISTACORE_HPP_

#include <vector>

#include "fistaParams.hpp"
#include "fdct3d.hpp"
#include "fdct3dinline.hpp"
#include "PARAMS.hpp"

class fistaCore{

	/**
	 * Public methods:
	 * Indeed, I encapsulated different types of migration in different public methods. The naming of those methods are self-explanatory.
	 **/

	public:
			void Reconstruct(fistaParams *fsp, PARAMS* params, fistaCore* fsc);
			void aprod(fistaParams* fsp, PARAMS* params, int operatorType);

	/**
	 * Access modifier of these methods are made private to prohibit the execution of these functions from anywhere except the object of this class
	 **/
	private:

			/* Read binary file. n1--> Fast direction (time), n2-> slower direction n3--> Slowest direction (space),
			 * data->data stored in this variable */
			int readBinFile(std::string fileName, int n1, int n2, int n3, std::vector<std::vector<std::vector<double> > > &data );

			int computeObsData(std::string fileName, int n1, int n2, int n3, std::vector<std::vector<std::vector<double> > > &data, std::vector<std::vector<int > > &sampleMat, std::vector<std::vector<std::vector<double> > > &obsData);

			void project(std::vector<std::vector<std::vector<double > > >x, double tau);

			inline int dotProduct(int n1, int n2, int n3, std::vector<std::vector<int> > &sampleMat, CpxNumTns &data );

			double norm(std::vector<std::vector<int > >&x, int normType);

			void wthresh(double thresh, std::vector< std::vector<CpxNumTns> > &x, std::vector<std::vector<CpxNumTns > >&y, std::vector<std::vector<std::vector<int> > > &cellStruct);

			inline void writeBinFile(std::string fileName, int n1, int n2, int n3, CpxNumTns &data);

};

#endif /* SPGLCORE_HPP_ */
