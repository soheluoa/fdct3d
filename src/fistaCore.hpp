/*
 * fistaCore.cpp
 *
 *  Created on: Sep 12, 2014
 *      Author: entropy
 */
/*  =================================================================================
 *  Name        : SeismicImaging.hpp
 *  Author      : Sohel Bhuiyan
 *  Version     : 0.1
 *  Purpose     : Base code of this software package. Contains all the basic methods entailed to generate image.
 *  			: This is a generalised migration code which can include different migration methods. It can accept
 *  			: any number of new method declaration.
 *  Date        : Feb 18, 2014
 *  Affiliation : University of Alberta, Physics department (SAIG)
 *  Email       : mbhuiyan@ualberta.ca
 * ===================================================================================
 */

#ifndef FISTACORE_HPP_
#define FISTACORE_HPP_

#include <vector>

#include "fistaParams.hpp"
#include "PARAMS.hpp"



class fistaCore{

	/**
	 * Public methods:
	 * Indeed, I encapsulated different types of migration in different public methods. The naming of those methods are self-explanatory.
	 **/

	public:
			void Execute(fistaParams *fsp, PARAMS* params, fistaCore* fsc);

	/**
	 * Access modifier of these methods are made private to prohibit the execution of these functions from anywhere except the object of this class
	 **/
	private:

			/* Read binary file. n1--> Fast direction (time), n2-> slower direction n3--> Slowest direction (space),
			 * data->data stored in this variable */
			int readBinFile(std::string fileName, int n1, int n2, int n3, std::vector<std::vector<std::vector<double> > > &data );

			int computeObsData(std::string fileName, int n1, int n2, int n3, std::vector<std::vector<std::vector<double> > > &data, std::vector<std::vector<int > > &sampleMat, std::vector<std::vector<std::vector<double> > > &obsData);

			/// x, and g are the curvelet/waveatom coefficients ///
			//void spgLineCurvy(double f, std::vector<std::vector<std::vector<double > > >x, std::vector<std::vector<std::vector<double > > >dx, double gtd, double fMax, std::vector<std::vector<std::vector<double > > >data);

			void project(std::vector<std::vector<std::vector<double > > >x, double tau);

			int spLine(double f, std::vector<std::vector<std::vector<double > > >x, std::vector<std::vector<std::vector<double > > >dx, double gtd, double fMax, std::vector<std::vector<std::vector<double> > > data);

			void aprod(fistaParams* fsp, PARAMS* params, int operatorType);

			double norm(std::vector<std::vector<int > >&x, int normType);

			void wthresh(int n1, int n2, int n3, vector<vector<vector<double> > > &data, double thresh, std::vector<std::vector<std::vector<double > > >&x);

			inline void writeBinFile(std::string fileName, int n1, int n2, int n3, std::vector<std::vector<std::vector<double> > > &data);

};

#endif /* SPGLCORE_HPP_ */
