/*
 * Spectra.h
 *
 *  Created on: 2011-02-23
 *      Author: weaver
 */

#ifndef SPECTRA_H_
#define SPECTRA_H_

#include "MultiAngleFLIC.h"
using namespace std;

class Spectra {

private:

public:

	double lambda [1000];
	double em [1000];

	Spectra() { } // private default constructor
	Spectra(string dye) {}

    int loadDye (string dye) {return 1;}

};

#endif /* SPECTRA_H_ */
