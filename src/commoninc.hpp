/*
    FDCT3D (Fast 3d Curvelet Transform)
	Written by Lexing Ying

	Updated by: Sohel Bhuiyan (2014)
				University of Alberta
	   	   	   	Upgraded from FFTW 2.1.5 to FFTW 3.3.3
	   	   	    Solved the memory leaking problem in FFTW 3.3.3
	   	   	    Plugged in with FISTA Algorithm
*/

#ifndef _FDCT3DINC_HPP_
#define _FDCT3DINC_HPP_

//STL stuff
#include <iostream>
#include <fstream>
#include <sstream>

#include <cfloat>
#include <cassert>
#include <cmath>
#include <string>
#include <string.h>
#include <complex>

#include <vector>
#include <set>
#include <map>
#include <deque>
#include <queue>
#include <utility>
#include <algorithm>
using namespace std;

//FFT stuff
#include "fftw3.h"

//typedef double double;
typedef complex<double> cpx;

//AUX functions
inline int pow2(int l) { assert(l>=0); return (1<<l); }


#endif
