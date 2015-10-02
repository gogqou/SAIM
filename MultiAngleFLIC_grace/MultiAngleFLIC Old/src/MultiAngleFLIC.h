/*
 * MultiAngleFLIC.h
 *
 *  Created on: 2010-09-04
 *      Author: matt
 */

#ifndef MULTIANGLEFLIC_H_
#define MULTIANGLEFLIC_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <complex>
#include <ctime>
#include <vector>
#include <list>
#include <iterator>
#include "mkl.h"
#include "mkl_rci.h"
#include "math.h"
#include <tiffio.h>
#include <limits>
#include <gsl/gsl_integration.h>
#include <Magick++.h>

using namespace std;

extern double DiI_Lam [];
extern double DiI_Em [];
extern double DiO_Lam [];
extern double DiO_Em [];

extern double P_Ex_Fromherz_MA(double theta2, double lambda, double dg);
extern double P_Em_Fromherz_Simple(double theta2, void* param);
extern double P_Em_theta(double lambda, void* param);
extern double P_Em_theta_lambda(double theta2Min, double theta2Max, double lambdaMin, double lambdaMax, double dg);
extern void generateFLICCurves(double* ZGAP, double* FLIC, double* DERIV, int nimages);
extern int nearestNeighbors(int* m, int* n, int i, int j, int nrows, int ncols);
extern int getIndex(int row, int col, int image, int nrows, int ncols);
extern int getCol(int index, int ncols);
extern int getRow(int index, int ncols);
extern double interpolateFLIC(double* curveX, double* curveY, double x_spacing, int sizeCurves, double x);
extern double interpolate(double* curveX, double* curveY, double x);
extern string IntToStr( int n );
extern int get_image_size(int first_image_number, int& nrows, int& ncols, string path, string image_name);
extern int importImages(int* images, int first_image_number, int nrows, int ncols, int nimages, string path, string image_name);
extern int importImage(int* Image, string file_path, int nrows, int ncols);
extern int arrayWriter(double* array, int elements, string file_path);
extern int arrayWriter(int* array, int elements, string file_path);
extern int readCSVtoArray(double* array, string file_path);
extern int nlsq_FLIC_local_HAB (double* Images, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double* A, double* B, double* R);
extern int nlsq_FLIC_global (double* Images, int* pixels, int num_pix, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double& A, double& B, double* R);
extern int nlsq_FLIC_Brute (double* Images, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double A, double B, double* R);
int nlsq_FLIC_BruteC (double* Images, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double A, double B, double* R, double cutoff, int* C);
extern int nlsq_FLIC_HAB_list (double* Images, int* pixels, int num_pix, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double& A, double& B, double* R);
extern void randPixels(int* pixels, int nrows, int ncols, int npixels);
extern void randPixelsC(int* pixels, double* Images, int nrows, int ncols, int npixels, double cutoff);
extern int loadImage(int* image, string file_path, int nrows, int ncols);
extern void getImageSize(string file_path, int& nrows, int& ncols);
extern void tester(double* x);
void testAFXN(void);


template <class T>
class csv_istream_iterator: public iterator<input_iterator_tag, T>
{
    istream * _input;
    char _delim;
    string _value;
public:
    csv_istream_iterator( char delim = ',' ): _input( 0 ), _delim( delim ) {}
    csv_istream_iterator( istream & in, char delim = ',' ): _input( &in ), _delim( delim ) { ++*this; }

    const T operator *() const {
        istringstream ss( _value );
        T value;
        ss >> value;
        return value;
    }

    istream & operator ++() {
        if( !( getline( *_input, _value, _delim ) ) )
        {
            _input = 0;
        }
        return *_input;
    }

    bool operator !=( const csv_istream_iterator & rhs ) const {
        return _input != rhs._input;
    }
};

template <>
const string csv_istream_iterator<string>::operator *() const {
    return _value;
}


#endif /* MULTIANGLEFLIC_H_ */
