/*
 * MultiAngleFLIC.h
 *
 *  Created on: 2010-09-04
 *      Author: matt
 */

#ifndef MULTIANGLEFLIC_H_
#define MULTIANGLEFLIC_H_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <complex>
#include <ctime>
#include <vector>
#include <list>
#include <iterator>
#include <algorithm>
#include <limits>
#include "mkl.h"
#include "mkl_spblas.h"
#include <tiffio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <Magick++.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "Spectra.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

double mean(double* x, int nvals);
double mean(int* x, int nvals);
extern void generate_image_seq_path(string& path, string image_path, string image_name, string file_ext, int image_num);
extern void generateAveragedRefSeries(string image_path, string base_name, int images_per_set, int nsets);
extern double P_Ex(double theta2, double lambda, double d2);
extern double P_Em(double theta2, void* param);
extern double P_Em_theta(double lambda, void* param);
extern double P_Em_theta_lambda(double theta2Min, double theta2Max, double lambdaMin, double lambdaMax, double d2);
extern double generateFLICCurves(double* ZGAP, double* FLIC, double* DERIV, int nimages);
extern int getNeighbors(int* m, int* n, int i, int j, int dist, int nrows, int ncols);
extern int getIndex(int row, int col, int image, int nrows, int ncols);
extern int getCol(int index, int ncols);
extern int getRow(int index, int ncols);
extern double interpolateFLIC(double* curveX, double* curveY, double x_spacing, int sizeCurves, double x);
extern double interpolate(double* curveX, double* curveY, double x);
extern string IntToStr( int n );
extern int StrToInt( const string& s );
extern double StrToDouble( const string& s );
extern int get_image_size(int first_image_number, int& nrows, int& ncols, string path, string image_name);
extern int importImages(int* images, int first_image_number, int nrows, int ncols, int nimages, string path, string image_name);
extern int importImage(int* Image, string file_path, int nrows, int ncols);
extern void writeImage(int* image, int width, int height, string output_file);
extern void writeImage(double* image, double multFactor, int width, int height, string output_file);
extern void writeImageJPEG(int* image, int width, int height, string output_file);
extern void writeImagePGM(int* image, int width, int height, string output_file);
extern int arrayWriter(double* array, int elements, string file_path);
extern int arrayWriter(int* array, int elements, string file_path);
extern int readCSVtoArray(double* array, string file_path);
extern int readINTCSVtoArray(int* array, string file_path);
extern void removeAllWhiteSpaces(string &str);
extern int localFLICfitter (double* Images, int nrows, int ncols, int nimages, double* ZGAP, double* FLIC, double* DERIV,
		double* H, double* A, double* B, double* L, double* R, int* F, string fit);
extern int localFLICfitterConstrained (double* Images, int pixIndex, int nrows, int ncols, int nimages, double* ZGAP, double* FLIC, double* DERIV,
		double* H, double* A, double* B, double* L, double* R, int* F, int constraint, double maxChange, string fit);
extern int nlsq_FLIC_HAB_global (double* Images, int* pixels, int num_pix, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double& A, double& B, double& R);
extern int nlsq_FLIC_HAB_pixel_MS (double* Images, int pixel, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double& H, double& A, double& B, double& R, double h_LB, double h_UB, double stepSize);
extern int nlsq_FLIC_HABL_pixel_MS (double* Images, int pixel, int nrows, int ncols, int nimages, double& H, double& A, double& B, double& L,
		double& R, double h_LB, double h_UB, double stepSize);
extern int nlsq_FLIC_HAB_pixel(double* Images, int pixel, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double& H, double& A, double& B, double& R, double h_LB, double h_UB);
extern int nlsq_FLIC_HAB_pixelNC (double* Images, int pixel, int nrows, int ncols, int nimages, double& H, double& A, double& B, double& R, double h_LB, double h_UB);
extern int nlsq_FLIC_HABL_pixelNC (double* Images, int pixel, int nrows, int ncols, int nimages, double& H, double& A, double& B, double& L, double& R, double h_LB, double h_UB);
extern int FLIC_fitter_Brute (double* Images, int *pixels, int npix, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double A, double B, double* R);
extern int FLIC_fitter_BruteT (double* Images, int* pixels, int npix, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double A, double B, double* R);
extern void randPixels(int* pixels, int nrows, int ncols, int npixels);
extern void randPixelsC(int* pixels, double* Images, int nrows, int ncols, int nimages, int ndesiredPix, double cutoff);
extern void allPixelsC(int* pixels, double* Images, int nrows, int ncols, int nimages, int& ndesiredPix, double cutoff);
extern void threshPixels(int* pixels, int* Thresh, int nrows, int ncols, int& nPix, int val);
extern int loadImage(int* image, string file_path, int nrows, int ncols);
extern void getImageSize(string file_path, int& nrows, int& ncols);
extern void tester(double* x);
void testAFXN(void);
extern int nlsq_FLIC_HAB_globalGradient (double* Images, int* fitPixels, int nfitPix, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double& A, double& B, double& R, double str);

extern int Surface(double* zgrid, double* x,double* y,double* z, MKL_INT npoints, double* xnodes,double* ynodes, MKL_INT nx,
		MKL_INT ny, double smoothness, char interp, char regularizer);

extern void Pex_Point_fcn (MKL_INT *m, MKL_INT *n, double *x, double *f);
extern void Pex_PointOx_fcn (MKL_INT *m, MKL_INT *n, double *x, double *f);
extern void Pex_CylinderHAB_fcn (MKL_INT *m, MKL_INT *n, double *x, double *f);
extern void Pex_SphereHAB_fcn (MKL_INT *m, MKL_INT *n, double *x, double *f);
extern void Pex_CylinderHABL_fcn (MKL_INT *m, MKL_INT *n, double *x, double *f);
extern void Pex_SphereHABL_fcn (MKL_INT *m, MKL_INT *n, double *x, double *f);

extern double Pex_SimpleOx(double thetaEx, void* param);
extern double Pex_SimpleOx_thetaEx(double thetaExMin, double thetaExMax, double theta2, double lambda, double dox);

extern void pixelSeriesGenerator(double* values, int nvals, string image_path, string image_name);


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
