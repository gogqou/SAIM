//============================================================================
// Name        : MultiAngleFLIC.cpp
// Author      : Matthew Paszek
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "MultiAngleFLIC.h"

using namespace std;

#define PI 3.14159265				//PI
#define DTR 2*PI/360				//Conversion factor from degrees to radians

const double n3 = 1.33;				//The refractive index of the ambient media / cytoplasm
const double n2 = 1.44;				//The refractive index of the cell membrane
const double ng = 1.33;				//The refractive index of the cell membrane / SiO2 gap
const double n1 = 1.46;				//The refractive index of SiO2
const double n0 = 4.00;				//The refractive index of Si
const double d2 = 4;				//The thickness of the cell membrane in units of nm
const double d1 = 500;				//The thickness of the SiO2 layer
const double thetaDye = PI/2;		//The dye excitation angle (usually pi/2 for DiI/DiO in cells)
const double lambdaEx = 488; //561;		//The wavelengDiI_Lamth of the microscope excitation laser
const double lambdaEmMin = 500; //580;		//The minimum wavelength of the microscope emission filter
const double lambdaEmMax = 550; //640;		//The maximum wavelength of the microscope emission filter
double* spectraEmX = DiO_Lam; //DiI_Lam;
double* spectraEmY = DiO_Em; //DiI_Em;
const double z_gap_max = 400;		//The maximum possible gap (glycocalyx) thickness (used for constructing the FLIC curve)
const double z_gap_spacing = 1;		//The desired sampling rate of points on the FLIC curve
const int FLIC_points = 			//The number of sampled points on the FLIC curve
		floor(z_gap_max/z_gap_spacing)+1;
const double derivative_step = 0.025;
const int num_images = 11;			//The number of images for each sample
double thetaLaser [] 		//The laser inclination angle (in radians) for each image
    = {0,5,10,15,20,25,30,35,40,45,50};
    //= {0,5,10,15,20,25,30,35,40,45,50,55,60}; //13
	//={0,6.66,13.33,20.00,26.66,33.33,40.00,46.66,53.33,60.00};  //10
	//={0,7.5,15,22.5,30,37.5,45,52.5,60};  //9
	//={0,8.57,17.14,25.71,34.28,42.85,51.42,59.99};  //8
	//={0,10,20,30,40,50,60};  //7
	//={0,12,24,36,48,60};  //6
	//={0,15,30,45,60};  //5
    //={0,20,40,60}; 	//4
    //={0,30,60};	//3


const double thetaObjMax = 67*DTR;  //Upper limit is equal to the critical angle from going from
									//membrane (n2) to medium (n3); Also must consider that the light with this
									//angle from the membrane can enter into the objective (dependent
									//on its NA) after refracting in the water

const int construct_FLIC_curves = 1;//A flag that indicates if the FLIC curve should be constructed (1) or read from file (!=1)
const int FLIC_image_start_num = 0;		//The first number of the image to be analyzed in a sequence
const int ref_image_start_num = 0;  //The first number of reference images
const int background_image_start_num = 0; //The first number of the background images
const int num_image_sets = 2;		//The number of sets of flic images per sample point
const double A_guess = MaxRGB*0.2;	//The initial guess of the scaling parameter for the height fitting algorithm
const double B_guess = MaxRGB*0.05;	//The initial guess of the background parameter for the height fitting algorithm
const double h_min = 10;			//The minimum possible height in the fitting algorithm
const double h_max = 390;			//The maximum possible height in the fitting algorithm
const double A_min = 30;			//The minimum possible scaling parameter in the fitting algorithm
const double A_max = MaxRGB*0.5;	//The maximum possible scaling parameter in the fitting algorithm
const double B_min = 0;				//The minimum possible background constant in the fitting algorithm
const double B_max = MaxRGB*0.1;	//The maximum possible background constant in the fitting algorithm
const double I_cutoff = MaxRGB*0.5;//The maximum pixel intensity.  Pixels with greater intensity are not processed
const int MS_step = 50;				//The sampling rate for initial height guesses in the multi-start method in the
									//fitting algorithm
const int nPixtoFit = 400;			//The number of randomly selected pixels that will be globally fit to calculate the best
									//fit A and B

/*
string FLIC_image_name ("100_4_02_01_");					//The base name of the main image files (should not include the image sequence #)
string background_image_name ("background_image_");		//The base name of the background image files
string ref_image_name ("reference_image_");			//The base name of the reference images
string ref_background_image_name ("ref_background_image");				//The base name of the background image for the reference images
string image_file_ext (".tif");						//The file extension on the image files.
string exp_name ("EXP100_4_02_01");									//The name of the experiment (used to name output files)
string image_path ("/home/matt/workspace/MultiAngleFLIC/images/");	//The path to the directory containing the image files
string FLIC_path ("/home/matt/workspace/MultiAngleFLIC/curves/");	//The path to the directory containing the FLIC curves
string output_path ("/home/matt/workspace/MultiAngleFLIC/output/");	//The path to the directory that the output files will be
string ZGAP_name ("ZGAP13.csv");
string FLIC_name ("FLIC13.csv");
string DERIV_name ("DERIV13.csv");
*/

string FLIC_image_name ("control");					//The base name of the main image files (should not include the image sequence #)
string background_image_name ("cont_background");		//The base name of the background image files
string ref_image_name ("reference");			//The base name of the reference images
string ref_background_image_name ("background");				//The base name of the background image for the reference images
string image_file_ext (".tif");						//The file extension on the image files.
string exp_name ("10ACont");									//The name of the experiment (used to name output files)
string image_path ("/home/matt/FLIC_for_analysis/processed/");	//The path to the directory containing the image files
string FLIC_path ("/home/matt/FLIC_for_analysis/processed/curves/");	//The path to the directory containing the FLIC curves
string output_path ("/home/matt/FLIC_for_analysis/processed/output/");	//The path to the directory that the output files will be
string ZGAP_name ("ZGAP_DiO_11.csv");
string FLIC_name ("FLIC_DiO_11.csv");
string DERIV_name ("DERIV_DiO_11.csv");


int main()
{
	//***********************************************The local variables**************************************************
	int i, j, k, counter, flag, num_rows, num_cols, index;
	double intensity;
	int *Images_temp;
	double *Images;
	double *H, *R;
	int *C;
	double A, B;
	double *ZGAP, *FLIC, *DERIV;
	int *globalPix;
	int nglobalPix;
	string file_path;

	//***********************************Construct the FLIC curves or read them from file**********************************

    //Allocate memory for the FLIC curves,
    ZGAP = new double[FLIC_points];
	FLIC = new double[FLIC_points*num_images];
	DERIV = new double[FLIC_points*num_images];

	//Convert the angles entered as degrees into radians;
	for (i=0; i<num_images; i++) thetaLaser[i] = thetaLaser[i]*DTR;

	//Generate the FLIC curves if the flag, construct_FLIC_Curves,  = 1, otherwise read in FLIC curves from file
	if (construct_FLIC_curves == 1) {

		//Generate the FLIC curves
		generateFLICCurves(ZGAP, FLIC, DERIV, num_images);  //construct the curves

		//Write the FLIC curves to file
		file_path = FLIC_path + ZGAP_name;
		arrayWriter(ZGAP, FLIC_points, file_path);
		file_path = FLIC_path + FLIC_name;
		arrayWriter(FLIC, FLIC_points*num_images, file_path);
		file_path = FLIC_path + DERIV_name;
		arrayWriter(DERIV, FLIC_points*num_images, file_path);

		return 0;
	}

	else {		//otherwise read in FLIC curves from file
		file_path = FLIC_path + ZGAP_name;
		if ( !readCSVtoArray(ZGAP, file_path) ) {
			cout << "\nFLIC curves (ZGAP.csv) could not be opened - Aborting program\n";
			return 1;
		}
		file_path = FLIC_path + FLIC_name;
		if ( !readCSVtoArray(FLIC, file_path) ) {
			cout << "\nFLIC curves (FLIC.csv) could not be opened - Aborting program\n";
			return 1;
		}
		file_path = FLIC_path + DERIV_name;
		if ( !readCSVtoArray(DERIV, file_path) ) {
			cout << "\nFLIC curves (DERIV.csv) could not be opened - Aborting program\n";
			return 1;
		}
	}

    //*****************************Open all the tiff images in the sequence and store****************************
	//*****************************their pixel values in row-major format in "images"****************************

	//Open the background images and record the average background for each of the main FLIC image
	if (background_image_start_num < 10)
		file_path = image_path + background_image_name + "00" + IntToStr( background_image_start_num ) + image_file_ext;
	else if (background_image_start_num < 100)
		file_path = image_path + background_image_name + "0" + IntToStr( background_image_start_num ) + image_file_ext;
	else file_path = image_path + background_image_name + IntToStr( background_image_start_num ) + image_file_ext;

	getImageSize(file_path, num_rows, num_cols);

	cout << "\nnum_rows: " << num_rows << " num_cols: " << num_cols << "\n";
	Images_temp = new int[ num_rows*num_cols*num_images*num_image_sets ];
	vector<double> background (num_images,0);
	intensity = 0;
	for (j = 0; j < num_image_sets; j++) {
		for (i = 0; i < num_images; i++) {
			index = j*num_images + i;

			//Set the image path
			if (index < 10 )
				file_path = image_path + background_image_name + "00" + IntToStr( background_image_start_num + index ) + image_file_ext;
			else if (index < 100)
				file_path = image_path + background_image_name + "0" + IntToStr( background_image_start_num + index ) + image_file_ext;
			else file_path = image_path + background_image_name + IntToStr( background_image_start_num + index ) + image_file_ext;

			//Load the image
			flag = loadImage(&Images_temp[index*num_rows*num_cols], file_path, num_rows, num_cols);  //import the images

			if (flag != 1) {
				cout << "\nAborting program - images in sequence are not all the same size";
				return 1;
			}

			//Record the average intensity in each of the background images
			intensity = 0;
			for (k = 0; k < num_rows*num_cols; k++) {
				intensity = intensity + Images_temp[k + index*num_rows*num_cols];
			}

			intensity = intensity/(num_rows*num_cols);
			background[i] = background[i] + intensity;
		}
	}

	for (i = 0; i < num_images; i++)
		background[i] = background[i]/num_image_sets;

	delete [] Images_temp;
	cout << "\nbackground[0]: " << background[0] << "\n";

	//Open the reference images and record the average intensity for each reference image
	if (ref_image_start_num < 10)
		file_path = image_path + ref_image_name + "00" + IntToStr( ref_image_start_num ) + image_file_ext;
	else if (ref_image_start_num < 100)
		file_path = image_path + ref_image_name + "0" + IntToStr( ref_image_start_num ) + image_file_ext;
	else
		file_path = image_path + ref_image_name + IntToStr( ref_image_start_num ) + image_file_ext;

	getImageSize(file_path, num_rows, num_cols);

	cout << "\nnum_rows: " << num_rows << " num_cols: " << num_cols << "\n";
	Images_temp = new int[ num_rows*num_cols*num_images*num_image_sets ];
	vector<double> reference (num_images,0);

	for (j = 0; j < num_image_sets; j++) {
		for (i = 0; i < num_images; i++) {

			index = j*num_images + i;

			//Set the image path
			if (ref_image_start_num + index < 10)
				file_path = image_path + ref_image_name + "00" + IntToStr( ref_image_start_num + index ) + image_file_ext;
			else if (ref_image_start_num + index < 100)
				file_path = image_path + ref_image_name + "0" + IntToStr( ref_image_start_num + index ) + image_file_ext;
			else
				file_path = image_path + ref_image_name + IntToStr( ref_image_start_num + index ) + image_file_ext;

			//Open the images
			flag = loadImage(&Images_temp[index*num_rows*num_cols], file_path, num_rows, num_cols);  //import the images

			if (flag != 1) {
				cout << "\nAborting program - images in sequence are not all the same size";
				return 1;
			}

			//Record the average intensity of the pixels in the reference image
			intensity = 0;
			for (k = 0; k < num_rows*num_cols; k++) {
				intensity = intensity + Images_temp[k + index*num_rows*num_cols];
			}

			intensity = intensity/(num_rows*num_cols);
			reference[i] = reference[i] + intensity;
		}
	}


	for (i = 0; i < num_images; i++)
		reference[i] = reference[i]/num_image_sets;

	delete [] Images_temp;
	cout << "reference[0]: " << reference[0] << "\n";

	//Open the background image for the reference slide and correct the reference intensities
	file_path = image_path + ref_background_image_name + image_file_ext;
	getImageSize(file_path, num_rows, num_cols);

	cout << "\nnum_rows: " << num_rows << " num_cols: " << num_cols << "\n";
	Images_temp = new int[ num_rows*num_cols ];
	double ref_background = 0;


	file_path = image_path + ref_background_image_name + image_file_ext;

	flag = loadImage(Images_temp, file_path, num_rows, num_cols);  //import the images

	if (flag != 1) {
		cout << "\nAborting program - images in sequence are not all the same size";
		return 1;
	}

	for (k = 0; k < num_rows*num_cols; k++) {
		ref_background = ref_background + Images_temp[k];
	}

	ref_background = ref_background/(num_rows*num_cols);

	for (k = 0; k<num_images; k++) {
		reference[k] = reference[k] - ref_background;
	}

	delete [] Images_temp;


	//Open the main FLIC image files
	if (FLIC_image_start_num < 10)
		file_path = image_path + FLIC_image_name + "00" + IntToStr( FLIC_image_start_num ) + image_file_ext;
	else if (FLIC_image_start_num < 100)
		file_path = image_path + FLIC_image_name + "0" + IntToStr( FLIC_image_start_num ) + image_file_ext;
	else
		file_path = image_path + FLIC_image_name + IntToStr( FLIC_image_start_num ) + image_file_ext;

	getImageSize(file_path, num_rows, num_cols);

	cout << "\nnum_rows: " << num_rows << " num_cols: " << num_cols << "\n";
	Images = new double[ num_rows*num_cols*num_images ];
	Images_temp = new int[ num_rows*num_cols*num_images*num_image_sets ];
	for (i = 0; i <  num_rows*num_cols*num_images; i++) Images[i] = 0;

	for (j = 0; j < num_image_sets; j++) {
		for (i = 0; i < num_images; i++) {

			index = j*num_images + i;

			if (FLIC_image_start_num + index < 10)
				file_path = image_path + FLIC_image_name + "00" + IntToStr( FLIC_image_start_num + index ) + image_file_ext;
			else if (FLIC_image_start_num + index < 100)
				file_path = image_path + FLIC_image_name + "0" + IntToStr( FLIC_image_start_num + index ) + image_file_ext;
			else file_path = image_path + FLIC_image_name + IntToStr( FLIC_image_start_num + index ) + image_file_ext;

			flag = loadImage(&Images_temp[index*num_rows*num_cols], file_path, num_rows, num_cols);  //import the images

			if (flag != 1) {
				cout << "\nAborting program - images in sequence are not all the same size";
				return 1;
			}
		}
	}

	//Load the average of the images in each image set into images
	for (i = 0; i < num_images; i++) {
		for (k = 0; k < num_rows*num_cols; k++) {
			Images[k + i*num_rows*num_cols];
		}
	}

	for (j = 0; j < num_image_sets; j++) {
		for (i = 0; i < num_images; i++) {
			for (k = 0; k < num_rows*num_cols; k++) {
				Images[k + i*num_rows*num_cols] = Images[k + i*num_rows*num_cols]
                          + Images_temp[k + i*num_rows*num_cols + j*num_rows*num_cols*num_images];
			}
		}
	}

	//Correct the intensities of the main FLIC images
	for (i = 0; i < num_images; i++) {
		for (j = 0; j < num_rows*num_cols; j++) {
			Images[i*num_rows*num_cols + j] =
					round( (Images[i*num_rows*num_cols + j]/num_image_sets - background[i])/reference[i] );
		}
	}


    //**********************Calculate the best-estimate fit membrane height at each pixel*************************

    //Allocate the matrix H for storing the best-estimate of the membrane height at each image pixel
    H = new double[ num_rows*num_cols ];
    R = new double[ num_rows*num_cols ];
    C = new int[ num_rows*num_cols ];

    //Get an estimate for A, B, and each pixel height by performing a local nlsq fit
    //for each pixel in a randomly selected list
    if (nPixtoFit > num_rows*num_cols) nglobalPix = num_rows*num_cols;
    else nglobalPix = nPixtoFit;

	globalPix = new int[nglobalPix];

	randPixelsC(globalPix, Images, num_rows, num_cols, nglobalPix, I_cutoff);

    nlsq_FLIC_HAB_list (Images, globalPix, nglobalPix, num_rows, num_cols, num_images, ZGAP,
    	FLIC, DERIV, H, A, B, R);

    //Perform a global optimization to get the best fit A and B for the randomly
    //selected pixels
    time_t start = time(NULL);
    nlsq_FLIC_global (Images, globalPix, nglobalPix, num_rows, num_cols, num_images, ZGAP,
        FLIC, DERIV, H, A, B, R);
    cout << "\ntime: " << time(NULL) - start << "\n";

    //Calculate the best fit height for every pixel in the image using the A and B
    //parameters acquired from the global fit.  Brute force is used to calculate the best-
    //fit heights
    nlsq_FLIC_BruteC (Images, num_rows, num_cols, num_images, ZGAP, FLIC, DERIV, H, A, B, R, I_cutoff, C);

    //Write the best fit heights to file
    file_path = output_path + exp_name + "_H.csv";
	arrayWriter(H, num_cols*num_rows, file_path);

	file_path = output_path + exp_name + "_R.csv";
	arrayWriter(R, num_cols*num_rows, file_path);

	file_path = output_path + exp_name + "_C.csv";
	arrayWriter(C, num_cols*num_rows, file_path);

	//output some info
	int max_h = 0;
	for (i = 0; i < num_rows*num_cols; i++) {
		if (H[i] > max_h) max_h = H[i];
	}
    cout << "\nA: " << A << " B: " << B << "max_h: " << max_h << "\n";

    //*********************************Delete dynamically allocated memory****************************************
    delete [] Images;
    delete [] ZGAP;
    delete [] FLIC;
    delete [] DERIV;
    delete [] H;
    delete [] R;
    delete [] C;
    delete [] globalPix;

    cout << "\nfinished!!!";
    return 0;
}



//The probability function for excitation or emission derived by Lambacher
//and Fromherz, 1996 and 2002 for a five layer system
double P_Ex_Fromherz_MA(double theta2, double lambda, double dg) {

	//Calculate each theta from theta2
	double thetaG = asin(n2*sin(theta2)/ng);
	double theta1 = asin(ng*sin(thetaG)/n1);
	double theta0 = asin(n1*sin(theta1)/n0);
	double theta3 = asin(n2*sin(theta2)/n3);

	//Calculate the phase shift
	double phi = (4*PI*n2/lambda)*d2*cos(theta2);

	//Calculate the pi's for calculating he fresnel coefficients
	double p0 = n0*cos(theta0);
	double p1 = n1*cos(theta1);
	double p2 = n2*cos(theta2);
	double pg = ng*cos(thetaG);

	//Calculate the qi's for calculating the fresnel coefficients
	double q0 = cos(theta0)/n0;
	double q1 = cos(theta1)/n1;
	double q2 = cos(theta2)/n2;
	double qg = cos(thetaG)/ng;

	//Calculate the fresnel coefficients
	double t32TE = 2*n3*cos(theta3)/(n3*cos(theta3) + n2*cos(theta2));
	double t32TM = 2*n3*cos(theta3)/(n2*cos(theta3) + n3*cos(theta2));
	double r23TE = (n2*cos(theta2) - n3*cos(theta3))/(n2*cos(theta2) + n3*cos(theta3));
	double r23TM = (n3*cos(theta2) - n2*cos(theta3))/(n3*cos(theta2) + n2*cos(theta3));

	//Calculate the Fresnesl coefficient for the layer system (gap, Si02, Si)
	double kg = 2*PI*ng/lambda;
	double k1 = 2*PI*n1/lambda;

	double lg = kg*dg*cos(thetaG);
	double l1 = k1*d1*cos(theta1);

	complex<double> i(0,1);

	complex<double> m11TE = cos(lg)*cos(l1) - (p1/pg)*sin(lg)*sin(l1);
	complex<double> m12TE = -(i/p1)*cos(lg)*sin(l1) - (i/pg)*sin(lg)*cos(l1);
	complex<double> m21TE = -i*pg*sin(lg)*cos(l1) - i*p1*cos(lg)*sin(l1);
	complex<double> m22TE = -(pg/p1)*sin(lg)*sin(l1) + cos(lg)*cos(l1);

	complex<double> m11TM = cos(lg)*cos(l1) - (q1/qg)*sin(lg)*sin(l1);
	complex<double> m12TM = (-i/q1)*cos(lg)*sin(l1) - (i/qg)*sin(lg)*cos(l1);
	complex<double> m21TM = -i*qg*sin(lg)*cos(l1) - i*q1*cos(lg)*sin(l1);
	complex<double> m22TM = -(qg/q1)*sin(lg)*sin(l1) + cos(lg)*cos(l1);

	complex<double> rTE = ((m11TE + m12TE*p0)*p2 - (m21TE + m22TE*p0)) /
			((m11TE + m12TE*p0)*p2 + (m21TE + m22TE*p0));

	complex<double> rTM = ((m11TM + m12TM*q0)*q2 - (m21TM + m22TM*q0)) /
			((m11TM + m12TM*q0)*q2 + (m21TM + m22TM*q0));


	//Calculate the probability of excitation or emmision

	complex<double> ifTE =   t32TE*(1. + rTE)/(1. - r23TE*rTE*exp(i*phi));

	complex<double> ifTM_p = t32TM*(1. - rTM)/(1. - r23TM*rTM*exp(i*phi));

	complex<double> ifTM_n = t32TM*(1. + rTM)/(1. - r23TM*rTM*exp(i*phi));

	double magSq_ifTE = abs(ifTE)*abs(ifTE);
	double magSq_ifTM_p = abs(ifTM_p)*abs(ifTM_p);
	double magSq_ifTM_n = abs(ifTM_n)*abs(ifTM_n);

	double P = sin(thetaDye)*sin(thetaDye)*magSq_ifTE +
			sin(thetaDye)*sin(thetaDye)*cos(theta2)*cos(theta2)*magSq_ifTM_p +
			2*cos(thetaDye)*cos(thetaDye)*sin(theta2)*sin(theta2)*magSq_ifTM_n;

	return P;

	//P could also be multiplied by an aperature function if necessary here

}

double P_Em_Fromherz_Simple(double theta2, void* param) {

	double* dbl_point;
	dbl_point = (double*) param;
	double lambda = dbl_point[0];
	double dg = dbl_point[1];

	//Calculate each theta from theta2
	double thetaG = asin(n2*sin(theta2)/ng);
	double theta1 = asin(ng*sin(thetaG)/n1);
	double theta0 = asin(n1*sin(theta1)/n0);
	double theta3 = asin(n2*sin(theta2)/n3);

	//Calculate the phase shift
	double phi = (4*PI*n2/lambda)*d2*cos(theta2);

	//Caculate the pi's for calculating he fresnel coefficients
	double p0 = n0*cos(theta0);
	double p1 = n1*cos(theta1);
	double p2 = n2*cos(theta2);
	double pg = ng*cos(thetaG);

	//Caclulate the qi's for calculating the fresnel coefficients
	double q0 = cos(theta0)/n0;
	double q1 = cos(theta1)/n1;
	double q2 = cos(theta2)/n2;
	double qg = cos(thetaG)/ng;

	//Calculate the fresnel coefficients
	double t23TE = 2*n2*cos(theta2)/(n2*cos(theta2) + n3*cos(theta3));
	double t23TM = 2*n2*cos(theta2)/(n3*cos(theta2) + n2*cos(theta3));
	double r23TE = (n2*cos(theta2) - n3*cos(theta3))/(n2*cos(theta2) + n3*cos(theta3));
	double r23TM = (n3*cos(theta2) - n2*cos(theta3))/(n3*cos(theta2) + n2*cos(theta3));

	//Calculate the Fresnesl coefficient for the layer system (gap, Si02, Si)
	double kg = 2*PI*ng/lambda;
	double k1 = 2*PI*n1/lambda;

	double lg = kg*dg*cos(thetaG);
	double l1 = k1*d1*cos(theta1);

	complex<double> i(0,1);

	complex<double> m11TE = cos(lg)*cos(l1) - (p1/pg)*sin(lg)*sin(l1);
	complex<double> m12TE = -(i/p1)*cos(lg)*sin(l1) - (i/pg)*sin(lg)*cos(l1);
	complex<double> m21TE = -i*pg*sin(lg)*cos(l1) - i*p1*cos(lg)*sin(l1);
	complex<double> m22TE = -(pg/p1)*sin(lg)*sin(l1) + cos(lg)*cos(l1);

	complex<double> m11TM = cos(lg)*cos(l1) - (q1/qg)*sin(lg)*sin(l1);
	complex<double> m12TM = (-i/q1)*cos(lg)*sin(l1) - (i/qg)*sin(lg)*cos(l1);
	complex<double> m21TM = -i*qg*sin(lg)*cos(l1) - i*q1*cos(lg)*sin(l1);
	complex<double> m22TM = -(qg/q1)*sin(lg)*sin(l1) + cos(lg)*cos(l1);

	complex<double> one = 1;
	complex<double> rTE = ((m11TE + m12TE*p0)*p2 - (m21TE + m22TE*p0)) /
			((m11TE + m12TE*p0)*p2 + (m21TE + m22TE*p0));

	complex<double> rTM = ((m11TM + m12TM*q0)*q2 - (m21TM + m22TM*q0)) /
			((m11TM + m12TM*q0)*q2 + (m21TM + m22TM*q0));


	//Calculate the probability of excitation or emission

	complex<double> ifTE =   t23TE*( 1. + rTE)/(1. - r23TE*rTE*exp(i*phi));

	complex<double> ifTM_p = t23TM*( 1. - rTM)/(1. - r23TM*rTM*exp(i*phi));

	complex<double> ifTM_n = t23TM*( 1. + rTM)/(1. - r23TM*rTM*exp(i*phi));

	double magSq_ifTE = abs(ifTE)*abs(ifTE);
	double magSq_ifTM_p = abs(ifTM_p)*abs(ifTM_p);
	double magSq_ifTM_n = abs(ifTM_n)*abs(ifTM_n);

	double U = sin(thetaDye)*sin(thetaDye)*magSq_ifTE +
			sin(thetaDye)*sin(thetaDye)*cos(theta2)*cos(theta2)*magSq_ifTM_p +
			2*cos(thetaDye)*cos(thetaDye)*sin(theta2)*sin(theta2)*magSq_ifTM_n;

	double Aout = n3*cos(theta3)/(n2*cos(theta2));

	double P = sin(theta2)*Aout*U;


	return P;
}

/*
double P_Em_Fromherz(double theta2, void* param) {

	double* dbl_point;
	dbl_point = (double*) param;
	double lambda = dbl_point[0];
	double dg = dbl_point[1];

	//Calculate each theta from theta2
	double thetaG = asin(n2*sin(theta2)/ng);
	double theta1 = asin(ng*sin(thetaG)/n1);
	double theta0 = asin(n1*sin(theta1)/n0);
	double theta3 = asin(n2*sin(theta2)/n3);

	//Calculate the phase shift
	double phi = (4*PI*n2/lambda)*d2*cos(theta2);

	//Caculate the pi's for calculating he fresnel coefficients
	double p0 = n0*cos(theta0);
	double p1 = n1*cos(theta1);
	double p2 = n2*cos(theta2);
	double pg = ng*cos(thetaG);

	//Caclulate the qi's for calculating the fresnel coefficients
	double q0 = cos(theta0)/n0;
	double q1 = cos(theta1)/n1;
	double q2 = cos(theta2)/n2;
	double qg = cos(thetaG)/ng;

	//Calculate the fresnel coefficients
	double t23TE = 2*n2*cos(theta2)/(n2*cos(theta2) + n3*cos(theta3));
	double t23TM = 2*n2*cos(theta2)/(n3*cos(theta2) + n2*cos(theta3));
	double r23TE = (n2*cos(theta2) - n3*cos(theta3))/(n2*cos(theta2) + n3*cos(theta3));
	double r23TM = (n3*cos(theta2) - n2*cos(theta3))/(n3*cos(theta2) + n2*cos(theta3));

	//Calculate the Fresnesl coefficient for the layer system (gap, Si02, Si)
	double kg = 2*PI*ng/lambda;
	double k1 = 2*PI*n1/lambda;

	double lg = kg*dg*cos(thetaG);
	double l1 = k1*d1*cos(theta1);

	complex<double> i(0,1);

	complex<double> m11TE = cos(lg)*cos(l1) - (p1/pg)*sin(lg)*sin(l1);
	complex<double> m12TE = -(i/p1)*cos(lg)*sin(l1) - (i/pg)*sin(lg)*cos(l1);
	complex<double> m21TE = -i*pg*sin(lg)*cos(l1) - i*p1*cos(lg)*sin(l1);
	complex<double> m22TE = -(pg/p1)*sin(lg)*sin(l1) + cos(lg)*cos(l1);

	complex<double> m11TM = cos(lg)*cos(l1) - (q1/qg)*sin(lg)*sin(l1);
	complex<double> m12TM = (-i/q1)*cos(lg)*sin(l1) - (i/qg)*sin(lg)*cos(l1);
	complex<double> m21TM = -i*qg*sin(lg)*cos(l1) - i*q1*cos(lg)*sin(l1);
	complex<double> m22TM = -(qg/q1)*sin(lg)*sin(l1) + cos(lg)*cos(l1);

	complex<double> one = 1;
	complex<double> rTE = ((m11TE + m12TE*p0)*p2 - (m21TE + m22TE*p0)) /
			((m11TE + m12TE*p0)*p2 + (m21TE + m22TE*p0));

	complex<double> rTM = ((m11TM + m12TM*q0)*q2 - (m21TM + m22TM*q0)) /
			((m11TM + m12TM*q0)*q2 + (m21TM + m22TM*q0));


	//Calculate the probability of excitation or emission

	complex<double> ifTE =   t23TE*( 1. + rTE)/(1. - r23TE*rTE*exp(i*phi));

	complex<double> ifTM_p = t23TM*( 1. - rTM)/(1. - r23TM*rTM*exp(i*phi));

	complex<double> ifTM_n = t23TM*( 1. + rTM)/(1. - r23TM*rTM*exp(i*phi));

	double magSq_ifTE = abs(ifTE)*abs(ifTE);
	double magSq_ifTM_p = abs(ifTM_p)*abs(ifTM_p);
	double magSq_ifTM_n = abs(ifTM_n)*abs(ifTM_n);

	double U = sin(thetaDye)*sin(thetaDye)*magSq_ifTE +
			sin(thetaDye)*sin(thetaDye)*cos(theta2)*cos(theta2)*magSq_ifTM_p +
			2*cos(thetaDye)*cos(thetaDye)*sin(theta2)*sin(theta2)*magSq_ifTM_n;

	double Aout = n3*cos(theta3)/(n2*cos(theta2));

	double P = sin(theta2)*Aout*U;


	return P;
}
*/


double P_Em_Fromherz(double theta2, void* param) {

	double* dbl_point;
	dbl_point = (double*) param;
	double lambda = dbl_point[0];
	double dg = dbl_point[1];

	//Calculate each theta from theta2
	double thetaG = asin(n2*sin(theta2)/ng);
	double theta1 = asin(ng*sin(thetaG)/n1);
	double theta0 = asin(n1*sin(theta1)/n0);
	double theta3 = asin(n2*sin(theta2)/n3);

	//Calculate the phase shift
	double phi = (4*PI*n2/lambda)*d2*cos(theta2);

	//Caculate the pi's for calculating he fresnel coefficients
	double p0 = n0*cos(theta0);
	double p1 = n1*cos(theta1);
	double p2 = n2*cos(theta2);
	double pg = ng*cos(thetaG);

	//Caclulate the qi's for calculating the fresnel coefficients
	double q0 = cos(theta0)/n0;
	double q1 = cos(theta1)/n1;
	double q2 = cos(theta2)/n2;
	double qg = cos(thetaG)/ng;

	//Calculate the fresnel coefficients
	double t23TE = 2*n2*cos(theta2)/(n2*cos(theta2) + n3*cos(theta3));
	double t23TM = 2*n2*cos(theta2)/(n3*cos(theta2) + n2*cos(theta3));
	double r23TE = (n2*cos(theta2) - n3*cos(theta3))/(n2*cos(theta2) + n3*cos(theta3));
	double r23TM = (n3*cos(theta2) - n2*cos(theta3))/(n3*cos(theta2) + n2*cos(theta3));

	//Calculate the Fresnesl coefficient for the layer system (gap, Si02, Si)
	double kg = 2*PI*ng/lambda;
	double k1 = 2*PI*n1/lambda;

	double lg = kg*dg*cos(thetaG);
	double l1 = k1*d1*cos(theta1);

	complex<double> i(0,1);

	complex<double> m11TE = cos(lg)*cos(l1) - (p1/pg)*sin(lg)*sin(l1);
	complex<double> m12TE = -(i/p1)*cos(lg)*sin(l1) - (i/pg)*sin(lg)*cos(l1);
	complex<double> m21TE = -i*pg*sin(lg)*cos(l1) - i*p1*cos(lg)*sin(l1);
	complex<double> m22TE = -(pg/p1)*sin(lg)*sin(l1) + cos(lg)*cos(l1);

	complex<double> m11TM = cos(lg)*cos(l1) - (q1/qg)*sin(lg)*sin(l1);
	complex<double> m12TM = (-i/q1)*cos(lg)*sin(l1) - (i/qg)*sin(lg)*cos(l1);
	complex<double> m21TM = -i*qg*sin(lg)*cos(l1) - i*q1*cos(lg)*sin(l1);
	complex<double> m22TM = -(qg/q1)*sin(lg)*sin(l1) + cos(lg)*cos(l1);

	complex<double> one = 1;
	complex<double> rTE = ((m11TE + m12TE*p0)*p2 - (m21TE + m22TE*p0)) /
			((m11TE + m12TE*p0)*p2 + (m21TE + m22TE*p0));

	complex<double> rTM = ((m11TM + m12TM*q0)*q2 - (m21TM + m22TM*q0)) /
			((m11TM + m12TM*q0)*q2 + (m21TM + m22TM*q0));


	//Calculate the probability of excitation or emission

	complex<double> ifTE =   t23TE*( 1. + rTE)/(1. - r23TE*rTE*exp(i*phi));

	complex<double> ifTM_p = t23TM*( 1. - rTM)/(1. - r23TM*rTM*exp(i*phi));

	complex<double> ifTM_n = t23TM*( 1. + rTM)/(1. - r23TM*rTM*exp(i*phi));

	double magSq_ifTE = abs(ifTE)*abs(ifTE);
	double magSq_ifTM_p = abs(ifTM_p)*abs(ifTM_p);
	double magSq_ifTM_n = abs(ifTM_n)*abs(ifTM_n);

	double U = sin(thetaDye)*sin(thetaDye)*magSq_ifTE +
			sin(thetaDye)*sin(thetaDye)*cos(theta2)*cos(theta2)*magSq_ifTM_p +
			2*cos(thetaDye)*cos(thetaDye)*sin(theta2)*sin(theta2)*magSq_ifTM_n;

	double Aout = n3*cos(theta3)/(n2*cos(theta2));

	double P = sin(theta2)*Aout*interpolate(spectraEmX, spectraEmY, lambda)*U;


	return P;
}


double P_Em_theta(double lambda, void* param) {


	double* dbl_point;
	dbl_point = (double*) param;
	double theta2Min = dbl_point[0];
	double theta2Max = dbl_point[1];
	double dg = dbl_point[2];

	double parameters [] = {lambda, dg};

	double result, error;

	gsl_integration_workspace * w
	         = gsl_integration_workspace_alloc (1000);


	gsl_function F;
	//F.function = &P_Em_Fromherz_Simple;
	F.function = &P_Em_Fromherz;
    F.params = &parameters;

	gsl_integration_qag(&F, theta2Min, theta2Max, 0, 1e-7, 1000, 3, w, &result, &error);

	gsl_integration_workspace_free (w);


	return result;
}

double P_Em_theta_lambda(double theta2Min, double theta2Max, double lambdaMin, double lambdaMax, double dg) {

	double param [] = {theta2Min, theta2Max, dg};
	double result, error;

	gsl_integration_workspace * w
			= gsl_integration_workspace_alloc (1000);

	gsl_function F;
	F.function = &P_Em_theta;
	F.params = &param;

	gsl_integration_qag(&F, lambdaMin, lambdaMax, 0, 1e-7, 1000, 3, w, &result, &error);

	gsl_integration_workspace_free (w);

	return result;
}

void generateFLICCurves(double* ZGAP, double* FLIC, double* DERIV, int nimages) {

	int i, j, counter;
	double Pex, max_I;
	double *Pem;

	//Allocate ZGAP, the z-heights of the membrane that the FLIC intensities will be calculated
	//for
	for (i = 0; i<FLIC_points; i++) {
		ZGAP[i] = i*z_gap_spacing;
	}

	//Calculate Pem for each height
	cout << "\nCalculating Pem for FLIC...\n";
	Pem = new double[FLIC_points];
	for (j=0; j < FLIC_points; j++) {
		Pem[j] = P_Em_theta_lambda(0, thetaObjMax, lambdaEmMin, lambdaEmMax, ZGAP[j]);
	}

	//Loop through all the input parameters and calculate Pex and Pex*Pem
	max_I = -1;
	counter = 0;
	for (i = 0; i < nimages; i++) {
		cout << "Generating FLIC curve: " << i << "\n";
	    for (j = 0; j < FLIC_points; j++) {
	        Pex = P_Ex_Fromherz_MA(thetaLaser[i], lambdaEx, ZGAP[j]);
	        FLIC[counter] = Pex*Pem[j];
	        if (FLIC[counter] > max_I) {
	        	max_I = FLIC[counter];
	        }
	        counter++;
	    }
	}

	//Calculate the derivatives of the FLIC curves
	cout << "\nCalculating Pem for DERIV...\n";
	for (j = 0; j < FLIC_points; j++) {
        Pem[j] = P_Em_theta_lambda(0, thetaObjMax, lambdaEmMin, lambdaEmMax, ZGAP[j]+derivative_step);
    }
    counter = 0;
	for (i = 0; i < nimages; i++) {
		cout << "Generating Deriv curve: " << i << "\n";
	    for (j = 0; j < FLIC_points; j++) {
	        Pex = P_Ex_Fromherz_MA(thetaLaser[i], lambdaEx, ZGAP[j]+derivative_step);
	        DERIV[counter] = (Pex*Pem[j] - FLIC[counter])/derivative_step;
	        counter++;
	    }
	}

	//Normalize the FLIC curves
	counter = 0;
	for (i = 0; i < nimages; i++) {
		    for (j = 0; j < FLIC_points; j++) {
				FLIC[counter] = FLIC[counter]/max_I;
				DERIV[counter] = DERIV[counter]/max_I;
				counter++;
		    }
	}
}


//fjac = df1/dx1, df2/dx1, df3/dx1, df4/dx1, df1/dx2, df2/dx2, df2/dx3, df2/dx4,...
/* 	nonlinear least square problem with boundary constraints

	Variables
	images	in:	TR solver handle
	heights out:     solution vector.  contains height of each pixel
	A       out:     solution for scaling factor A
	B       out:     solution for background parameter B
*/
int nlsq_FLIC_local_HAB (double* Images, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double* A, double* B, double* R)
{
	//***********************************Variables, etc*******************************************
	MKL_INT n = 3;			//number of function variables (height of pixel, A, B)
	MKL_INT m = nimages;	//dimension of function value (number of equations)
	double	eps[6];			// precisions for stop-criteria (see manual for more details
	double	*x;				// solution vector. contains values x for f(x)
	MKL_INT	iter1 = 1000;	//maximum number of iterations
	MKL_INT	iter2 = 100;	//maximum number of iterations of calculation of trial-step
	double	rs = 0.0;		//initial step bound
	MKL_INT RCI_Request;	//reverse communication interface parameter
	MKL_INT successful;		//controls of rci cycle
	double *fvec;			//function (f(x)) value vector
	double *fobs;			//The observed value for each function
	double *fjac;			//jacobi matrix
	double	*LW, *UP;		//lower and upper bounds
	MKL_INT iter;			//number of iterations
	MKL_INT st_cr;			//number of stop-criterion
	double r1, r2;			//initial and final residuals
	_TRNSPBC_HANDLE_t handle;	//TR solver handle
	int *NN_row, *NN_col;	//Storage of (row, col) of nearest neighbor pixels
	int num_NN;				//The number of neighbors for a given pixel
	MKL_INT row, col, i, j, k;		//counters
	double h_best, A_best, B_best;	//The best multi-start solution for h, A, and B
	double r_min, r1_min;			//Holds the lowest residual from the multi-start solutions
	double MS_start;

	//memory allocation
	x = new double[n];
	fvec = new double[m];
	fobs = new double[m];
	fjac = new double[m*n];
	LW = new double[n];
	UP = new double[n];
	NN_row = new int[9];
	NN_col = new int[9];

	//set precisions for stop-criteria
	for (i = 0; i < 6; i++)
	{
		eps [i] = 0.0000001;
	}

	//set bounds
	LW [0] = h_min;
	LW [1] = A_min;
	LW [2] = B_min;
	UP [0] = h_max;
	UP [1] = A_max;
	UP [2] = B_max;

	//Loop through all the pixels and calculate the best estimate of the height for each pixel
	for (row = 0; row < nrows; row++) {
		for (col = 0; col < ncols; col++) {

			//Get the nearest neighbor pixels
			num_NN = nearestNeighbors(NN_row, NN_col, row, col, nrows, ncols);

			//Compute fobs
			for (i = 0; i < num_images; i++) {
				for (j = 0; j < num_NN; j++) {
					fobs[i] = fobs[i] + Images[getIndex(NN_row[j], NN_col[j], i, nrows, ncols)];
				}
				fobs[i] = fobs[i]/num_NN;
			}

			//Start a multi-start approach to find the best height
			r_min = numeric_limits<double>::infinity();
			r1_min = numeric_limits<double>::infinity();

			for (MS_start = h_min; MS_start <= h_max; MS_start = MS_start + MS_step) {
				//Set the initial guess
				x[0] = MS_start;
				x[1] = A_guess;
				x[2] = B_guess;

				//set initial values
				for (i = 0; i < m; i++)
					fvec [i] = 0.0;
				for (i = 0; i < m*n; i++)
					fjac [i] = 0.0;


				//***********************Solve for the best height, A, and B********************************

				/* initialize solver (allocate mamory, set initial values)
					handle	in/out:	TR solver handle
					n       in:     number of function variables
					m       in:     dimension of function value
					x       in:     solution vector. contains values x for f(x)
					LW		in:		lower bound
					UP		in:		upper bound
					eps     in:     precisions for stop-criteria
					iter1   in:     maximum number of iterations
					iter2   in:     maximum number of iterations of calculation of trial-step
					rs      in:     initial step bound */
				if (dtrnlspbc_init (&handle, &n, &m, x, LW, UP, eps, &iter1, &iter2, &rs) !=
					TR_SUCCESS)
				{
					// Exit if unsuccessful
					cout << "| error in dtrnlspbc_init\n";
					MKL_FreeBuffers();
					return 1;
				}

				//set initial rci cycle variables
				RCI_Request = 0;
				successful = 0;

				//rci cycle
				while (successful == 0)
				{
					/* call tr solver
						handle		in/out:	tr solver handle
						fvec		in:     vector
						fjac		in:     jacobi matrix
						RCI_request in/out:	return number which denote next step for performing */
					if (dtrnlspbc_solve (&handle, fvec, fjac, &RCI_Request) != TR_SUCCESS)
					{
						// Exit if unsuccessful
						cout << "| error in dtrnlspbc_solve\n";
						MKL_FreeBuffers();
						return 1;
					}

					// according with rci_request value we do next step
					if (RCI_Request == -1 ||
						RCI_Request == -2 ||
						RCI_Request == -3 ||
						RCI_Request == -4 ||
						RCI_Request == -5 ||
						RCI_Request == -6)
						// exit rci cycle
						successful = 1;

					if (RCI_Request == 1)
					{
						// recalculate function value
						for (i = 0; i < m; i++) {
							fvec[i] = x[1]*interpolateFLIC(ZGAP, &FLIC[i*FLIC_points], z_gap_spacing, FLIC_points, x[0]) + x[2] - fobs[i];
						}
					}
					if (RCI_Request == 2)
					{
						// compute jacobi matrix
						//fjac = df1/dx1, df2/dx1, df3/dx1, df4/dx1,..., df1/dx2, df2/dx2, df2/dx2, df2/dx2,...
						for (i = 0; i < m; i++) {
							fjac[i] = x[1]*interpolateFLIC(ZGAP, &DERIV[i*FLIC_points], z_gap_spacing, FLIC_points, x[0]);
							fjac[m + i] = interpolateFLIC(ZGAP, &FLIC[i*FLIC_points], z_gap_spacing, FLIC_points, x[0]);
							fjac[2*m + i] = 1;
						}
					}
				}

				/* get solution statuses
					handle            in:	TR solver handle
					iter              out:	number of iterations
					st_cr             out:	number of stop criterion
					r1                out:	initial residuals
					r2                out:	final residuals */
				if (dtrnlspbc_get (&handle, &iter, &st_cr, &r1, &r2) != TR_SUCCESS)
				{
					//exit if unsuccessful
					cout << "| error in dtrnlspbc_get\n";
					MKL_FreeBuffers();
					return 1;
				}

				// free handle memory
				if (dtrnlspbc_delete (&handle) != TR_SUCCESS)
				{
					//exit if unsuccessful
					cout << "| error in dtrnlspbc_delete\n";
					MKL_FreeBuffers();
					return 1;
				}

				// Release internal MKL memory
				MKL_FreeBuffers();

				//******************End of Solve for the best height, A, and B****************************

				//Check to see if the solution is better than previous multi-start solutions
				if (r2 < r_min) {
					h_best = x[0];
					A_best = x[1];
					B_best = x[2];
					r_min = r2;
				}
				if (r1 < r1_min) {
					r1_min = r1;
				}
			}
			H[ getIndex(row, col, 0, nrows, ncols) ] = h_best;
			A[ getIndex(row, col, 0, nrows, ncols) ] = A_best;
			B[ getIndex(row, col, 0, nrows, ncols) ] = B_best;
			R[ getIndex(row, col, 0, nrows, ncols) ] = r_min;
		}
	}


	/* free allocated memory */
	delete [] x;
	delete [] fvec;
	delete [] fjac;
	delete [] LW;
	delete [] UP;

	return 1;

}

int nlsq_FLIC_HAB_list (double* Images, int* pixels, int npix, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double& A, double& B, double* R)
{
	//***********************************Variables, etc*******************************************
	MKL_INT n = 3;			//number of function variables (height of pixel, A, B)
	MKL_INT m = nimages;	//dimension of function value (number of equations)
	double	eps[6];			// precisions for stop-criteria (see manual for more details
	double	*x;				// solution vector. contains values x for f(x)
	MKL_INT	iter1 = 1000;	//maximum number of iterations
	MKL_INT	iter2 = 100;	//maximum number of iterations of calculation of trial-step
	double	rs = 0.0;		//initial step bound
	MKL_INT RCI_Request;	//reverse communication interface parameter
	MKL_INT successful;		//controls of rci cycle
	double *fvec;			//function (f(x)) value vector
	double *fobs;			//The observed value for each function
	double *fjac;			//jacobi matrix
	double	*LW, *UP;		//lower and upper bounds
	MKL_INT iter;			//number of iterations
	MKL_INT st_cr;			//number of stop-criterion
	double r1, r2;			//initial and final residuals
	_TRNSPBC_HANDLE_t handle;	//TR solver handle
	int *NN_row, *NN_col;	//Storage of (row, col) of nearest neighbor pixels
	int num_NN;				//The number of neighbors for a given pixel
	MKL_INT row, col, i, j, k;		//counters
	double h_best, A_best, B_best;	//The best multi-start solution for h, A, and B
	double r_min, r1_min;			//Holds the lowest residual from the multi-start solutions
	double MS_start;

	//memory allocation
	x = new double[n];
	fvec = new double[m];
	fobs = new double[m];
	fjac = new double[m*n];
	LW = new double[n];
	UP = new double[n];
	NN_row = new int[9];
	NN_col = new int[9];

	//set precisions for stop-criteria
	for (i = 0; i < 6; i++)
	{
		eps [i] = 0.0000001;
	}

	//set bounds
	LW [0] = h_min;
	LW [1] = A_min;
	LW [2] = B_min;
	UP [0] = h_max;
	UP [1] = A_max;
	UP [2] = B_max;

	//Initialize A and B
	A = 0;
	B = 0;

	//Loop through all the pixels and calculate the best estimate of the height for each pixel
	for (k = 0; k < npix; k++) {

		row = getRow(pixels[k], ncols);
		col = getCol(pixels[k], ncols);

		//Get the nearest neighbor pixels
		num_NN = nearestNeighbors(NN_row, NN_col, row, col, nrows, ncols);

		//Compute fobs
		for (i = 0; i < num_images; i++) {
			for (j = 0; j < num_NN; j++) {
				fobs[i] = fobs[i] + Images[getIndex(NN_row[j], NN_col[j], i, nrows, ncols)];
			}
			fobs[i] = fobs[i]/num_NN;
		}

		//Start a multi-start approach to find the best height
		r_min = numeric_limits<double>::infinity();

		for (MS_start = h_min; MS_start <= h_max; MS_start = MS_start + MS_step) {
			//Set the initial guess
			x[0] = MS_start;
			x[1] = A_guess;
			x[2] = B_guess;

			//set initial values
			for (i = 0; i < m; i++)
				fvec [i] = 0.0;
			for (i = 0; i < m*n; i++)
				fjac [i] = 0.0;


			//***********************Solve for the best height, A, and B********************************

			/* initialize solver (allocate mamory, set initial values)
				handle	in/out:	TR solver handle
				n       in:     number of function variables
				m       in:     dimension of function value
				x       in:     solution vector. contains values x for f(x)
				LW		in:		lower bound
				UP		in:		upper bound
				eps     in:     precisions for stop-criteria
				iter1   in:     maximum number of iterations
				iter2   in:     maximum number of iterations of calculation of trial-step
				rs      in:     initial step bound */
			if (dtrnlspbc_init (&handle, &n, &m, x, LW, UP, eps, &iter1, &iter2, &rs) !=
				TR_SUCCESS)
			{
				// Exit if unsuccessful
				cout << "| error in dtrnlspbc_init\n";
				MKL_FreeBuffers();
				return 1;
			}

			//set initial rci cycle variables
			RCI_Request = 0;
			successful = 0;

			//rci cycle
			while (successful == 0)
			{
				/* call tr solver
					handle		in/out:	tr solver handle
					fvec		in:     vector
					fjac		in:     jacobi matrix
					RCI_request in/out:	return number which denote next step for performing */
				if (dtrnlspbc_solve (&handle, fvec, fjac, &RCI_Request) != TR_SUCCESS)
				{
					// Exit if unsuccessful
					cout << "| error in dtrnlspbc_solve\n";
					MKL_FreeBuffers();
					return 1;
				}

				// according with rci_request value we do next step
				if (RCI_Request == -1 ||
					RCI_Request == -2 ||
					RCI_Request == -3 ||
					RCI_Request == -4 ||
					RCI_Request == -5 ||
					RCI_Request == -6)
					// exit rci cycle
					successful = 1;

				if (RCI_Request == 1)
				{
					// recalculate function value
					for (i = 0; i < m; i++) {
						fvec[i] = x[1]*interpolateFLIC(ZGAP, &FLIC[i*FLIC_points], z_gap_spacing, FLIC_points, x[0]) + x[2] - fobs[i];
					}
				}
				if (RCI_Request == 2)
				{
					// compute jacobi matrix
					//fjac = df1/dx1, df2/dx1, df3/dx1, df4/dx1,..., df1/dx2, df2/dx2, df2/dx2, df2/dx2,...
					for (i = 0; i < m; i++) {
						fjac[i] = x[1]*interpolateFLIC(ZGAP, &DERIV[i*FLIC_points], z_gap_spacing, FLIC_points, x[0]);
						fjac[m + i] = interpolateFLIC(ZGAP, &FLIC[i*FLIC_points], z_gap_spacing, FLIC_points, x[0]);
						fjac[2*m + i] = 1;
					}
				}
			}

			/* get solution statuses
				handle            in:	TR solver handle
				iter              out:	number of iterations
				st_cr             out:	number of stop criterion
				r1                out:	initial residuals
				r2                out:	final residuals */
			if (dtrnlspbc_get (&handle, &iter, &st_cr, &r1, &r2) != TR_SUCCESS)
			{
				//exit if unsuccessful
				cout << "| error in dtrnlspbc_get\n";
				MKL_FreeBuffers();
				return 1;
			}

			// free handle memory
			if (dtrnlspbc_delete (&handle) != TR_SUCCESS)
			{
				//exit if unsuccessful
				cout << "| error in dtrnlspbc_delete\n";
				MKL_FreeBuffers();
				return 1;
			}

			// Release internal MKL memory
			MKL_FreeBuffers();

			//******************End of Solve for the best height, A, and B****************************

			//Check to see if the solution is better than previous multi-start solutions
			if (r2 < r_min) {
				h_best = x[0];
				A_best = x[1];
				B_best = x[2];
				r_min = r2;
			}
		}
		H[k] = h_best;
		R[k] = r_min;
		A = A + A_best;
		B = B + B_best;
		cout << "\nhbest: " << h_best << " A_best: " << A_best << " B_best: " << B_best << " r2: " << r_min << " iter " << iter;
	}

	A = A/npix;
	B = B/npix;

	/* free allocated memory */
	delete [] x;
	delete [] fvec;
	delete [] fjac;
	delete [] LW;
	delete [] UP;

	return 1;

}


int nlsq_FLIC_Brute (double* Images, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double A, double B, double* R)
{
	//***********************************Variables, etc*******************************************
	int i, j, row, col;		//iterators
	int num_NN;				//the number of nearest neighbors for a given pixel
	int NN_row [9];			//the row numbers of the nearest neighbors
	int NN_col [9];			//the column numbers of the nearest neighbors
	int *fobs =
		new int[nimages];	//the observed data
	double r;				//the residual
	double r_min;			//the minimum residual found
	double h_best;			//the height for a pixel with the lowest residual

	//Loop through all the pixels and calculate the best estimate of the height for each pixel
	for (row = 0; row < nrows; row++) {
		for (col = 0; col < ncols; col++) {

			//Get the nearest neighbor pixels
			num_NN = nearestNeighbors(NN_row, NN_col, row, col, nrows, ncols);

			//Compute fobs
			for (i = 0; i < num_images; i++) {
				for (j = 0; j < num_NN; j++) {
					fobs[i] = fobs[i] + Images[getIndex(NN_row[j], NN_col[j], i, nrows, ncols)];
				}
				fobs[i] = fobs[i]/num_NN;
			}

			//Find the best height for each pixel with brute force
			r_min = numeric_limits<double>::infinity();

			for (j = 0; ZGAP[j] < h_max; j++) {

				//Compute the residual
				r = 0;
				for (i = 0; i < nimages; i++) {
					r = r + (A*FLIC[i*FLIC_points + j] + B - fobs[i])*(A*FLIC[i*FLIC_points + j] + B - fobs[i]);
				}

				//Check to see if the solution is better than previous multi-start solutions
				if (r < r_min) {
					h_best = ZGAP[j];
					r_min = r;
				}
			}

			H[ getIndex(row, col, 0, nrows, ncols) ] = h_best;
			R[ getIndex(row, col, 0, nrows, ncols) ] = r_min;
		}
	}

	delete [] fobs;
	return 1;

}

int nlsq_FLIC_BruteC (double* Images, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double A, double B, double* R, double cutoff, int* C)
{
	//***********************************Variables, etc*******************************************
	int i, j, row, col;		//iterators
	int num_NN;				//the number of nearest neighbors for a given pixel
	int NN_row [9];			//the row numbers of the nearest neighbors
	int NN_col [9];			//the column numbers of the nearest neighbors
	int *fobs =
		new int[nimages];	//the observed data
	double r;				//the residual
	double r_min;			//the minimum residual found
	double h_best;			//the height for a pixel with the lowest residual
	int process;

	//Loop through all the pixels and calculate the best estimate of the height for each pixel
	for (row = 0; row < nrows; row++) {
		for (col = 0; col < ncols; col++) {

			//Get the nearest neighbor pixels
			num_NN = nearestNeighbors(NN_row, NN_col, row, col, nrows, ncols);

			//Make sure none of the pixels is above the cutoff threshold
			process = 1;
			for (j = 0; j < num_NN; j++) {
				if ( Images[getIndex(NN_row[j], NN_col[j], 0, nrows, ncols)] > I_cutoff) {
					C[ getIndex(row, col, 0, nrows, ncols) ] = 0;
					H[ getIndex(row, col, 0, nrows, ncols) ] = -1;
					R[ getIndex(row, col, 0, nrows, ncols) ] = -1;
					process = 0;
					break;
				}
			}

			if (process) {
				//Compute fobs
				for (i = 0; i < num_images; i++) {
					for (j = 0; j < num_NN; j++) {
						fobs[i] = fobs[i] + Images[getIndex(NN_row[j], NN_col[j], i, nrows, ncols)];
					}
					fobs[i] = fobs[i]/num_NN;
				}

				//Find the best height for each pixel with brute force
				r_min = numeric_limits<double>::infinity();

				for (j = 0; ZGAP[j] < h_max; j++) {

					//Compute the residual
					r = 0;
					for (i = 0; i < nimages; i++) {
						r = r + (A*FLIC[i*FLIC_points + j] + B - fobs[i])*(A*FLIC[i*FLIC_points + j] + B - fobs[i]);
					}

					//Check to see if the solution is better than previous multi-start solutions
					if (r < r_min) {
						h_best = ZGAP[j];
						r_min = r;
					}
				}

				H[ getIndex(row, col, 0, nrows, ncols) ] = h_best;
				R[ getIndex(row, col, 0, nrows, ncols) ] = r_min;
				C[ getIndex(row, col, 0, nrows, ncols) ] = 1;
			}
		}
	}

	delete [] fobs;
	return 1;

}



int nlsq_FLIC_global (double* Images, int* pixels, int npix, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double& A, double& B, double* R)
{
	//***********************************Variables, etc*******************************************
	MKL_INT n = npix+2;			//number of function variables (height of pixel, A, B)
	MKL_INT m = npix*nimages;	//dimension of function value (number of equations)
	double	eps[6];			// precisions for stop-criteria (see manual for more details
	double	*x;				// solution vector. contains values x for f(x)
	MKL_INT	iter1 = 1000;	//maximum number of iterations
	MKL_INT	iter2 = 100;	//maximum number of iterations of calculation of trial-step
	double	rs = 0.0;		//initial step bound
	MKL_INT RCI_Request;	//reverse communication interface parameter
	MKL_INT successful;		//controls of rci cycle
	double *fvec;			//function (f(x)) value vector
	double *fobs;			//The observed value for each function
	double *fjac;			//jacobi matrix
	double	*LW, *UP;		//lower and upper bounds
	MKL_INT iter;			//number of iterations
	MKL_INT st_cr;			//number of stop-criterion
	double r1, r2;			//initial and final residuals
	_TRNSPBC_HANDLE_t handle;	//TR solver handle
	MKL_INT row, col, i, j, k;		//counters
	double h_best, A_best, B_best;	//The best multi-start solution for h, A, and B
	int counter;


	//memory allocation
	x = new double[n];
	fvec = new double[m];
	fobs = new double[m];
	fjac = new double[m*n];
	LW = new double[n];
	UP = new double[n];

	//set precisions for stop-criteria
	for (i = 0; i < 6; i++)
	{
		eps [i] = 0.0000001;
	}

	//set bounds
	for (i = 0; i < npix; i++) {
		LW [i] = h_min;
		UP [i] = h_max;
	}
	LW [npix] 	= A_min;
	LW [npix+1] = B_min;
	UP [npix] 	= A_max;
	UP [npix+1] = B_max;

	//Set the initial guess
	for (i = 0; i < npix; i++) x[i] = H[i];
	x[npix]   = A;
	x[npix+1] = B;

	cout << "AB: " << A << " " << B;

	//set initial values
	for (i = 0; i < m; i++)
		fvec [i] = 0.0;
	for (i = 0; i < m*n; i++)
		fjac [i] = 0.0;

	//Set fobs
	counter = 0;
	for (i = 0; i < npix; i++) {
		row = getRow(pixels[i], ncols);
		col = getCol(pixels[i], ncols);
		for (j = 0; j < nimages; j++) {
			fobs[counter] = Images[getIndex(row, col, j, nrows, ncols)];
			counter++;
		}
	}



	//***********************Solve for the best height, A, and B********************************

	/* initialize solver (allocate mamory, set initial values)
		handle	in/out:	TR solver handle
		n       in:     number of function variables
		m       in:     dimension of function value
		x       in:     solution vector. contains values x for f(x)
		LW		in:		lower bound
		UP		in:		upper bound
		eps     in:     precisions for stop-criteria
		iter1   in:     maximum number of iterations
		iter2   in:     maximum number of iterations of calculation of trial-step
		rs      in:     initial step bound */
	if (dtrnlspbc_init (&handle, &n, &m, x, LW, UP, eps, &iter1, &iter2, &rs) != TR_SUCCESS)
	{
		// Exit if unsuccessful
		cout << "| error in dtrnlspbc_init\n";
		MKL_FreeBuffers();
		return 1;
	}

	//set initial rci cycle variables
	RCI_Request = 0;
	successful = 0;

	//rci cycle
	while (successful == 0)
	{
		/* call tr solver
			handle		in/out:	tr solver handle
			fvec		in:     vector
			fjac		in:     jacobi matrix
			RCI_request in/out:	return number which denote next step for performing */
		if (dtrnlspbc_solve (&handle, fvec, fjac, &RCI_Request) != TR_SUCCESS)
		{
			// Exit if unsuccessful
			cout << "| error in dtrnlspbc_solve\n";
			MKL_FreeBuffers();
			return 1;
		}

		// according with rci_request value we do next step
		if (RCI_Request == -1 ||
			RCI_Request == -2 ||
			RCI_Request == -3 ||
			RCI_Request == -4 ||
			RCI_Request == -5 ||
			RCI_Request == -6)
			// exit rci cycle
			successful = 1;

		if (RCI_Request == 1)
		{
			// recalculate function value
			for (i = 0; i < npix; i++) {
				for (j = 0; j < nimages; j++) {
					fvec[i*nimages + j] = x[npix]*interpolateFLIC(ZGAP, &FLIC[j*FLIC_points], z_gap_spacing, FLIC_points, x[i]) + x[npix+1] - fobs[i*nimages + j];
				}
			}
		}
		if (RCI_Request == 2)
		{
			// compute jacobi matrix
			//fjac = df1/dx1, df2/dx1, df3/dx1, df4/dx1,..., df1/dx2, df2/dx2, df2/dx2, df2/dx2,...
			for (i = 0; i < npix; i++) {
				for (j = 0; j < nimages; j++) {
					fjac[i*m + i*nimages + j] = x[npix]*interpolateFLIC(ZGAP, &DERIV[j*FLIC_points], z_gap_spacing, FLIC_points, x[i]);
					fjac[npix*m + i*nimages + j] = interpolateFLIC(ZGAP, &FLIC[j*FLIC_points], z_gap_spacing, FLIC_points, x[i]);
					fjac[(npix+1)*m + i*nimages + j] = 1;
					//fjac[(npix+2)*m + i*nimages + j] = ...;  if dye dipole angle is also fit
				}
			}
		}
	}

	/* get solution statuses
		handle            in:	TR solver handle
		iter              out:	number of iterations
		st_cr             out:	number of stop criterion
		r1                out:	initial residuals
		r2                out:	final residuals */
	if (dtrnlspbc_get (&handle, &iter, &st_cr, &r1, &r2) != TR_SUCCESS)
	{
		//exit if unsuccessful
		cout << "| error in dtrnlspbc_get\n";
		MKL_FreeBuffers();
		return 1;
	}

	// free handle memory
	if (dtrnlspbc_delete (&handle) != TR_SUCCESS)
	{
		//exit if unsuccessful
		cout << "| error in dtrnlspbc_delete\n";
		MKL_FreeBuffers();
		return 1;
	}

	// Release internal MKL memory
	MKL_FreeBuffers();

	//******************End of Solve for the best height, A, and B****************************

	//Copy the best-fit variables to the output variables
	for (i = 0; i < npix; i++) {
		H[i] = x[i];
	}
	A = x[npix];
	B = x[npix+1];


	/* free allocated memory */
	delete [] x;
	delete [] fvec;
	delete [] fjac;
	delete [] LW;
	delete [] UP;

	return 1;

}


int nlsq_FLIC_globalA (double* Images, int* pixels, int npix, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double* Ai, double& A, double& B, double* R)
{
	//***********************************Variables, etc*******************************************
	MKL_INT n = npix+nimages+2;			//number of function variables (height of pixel, A, B)
	MKL_INT m = npix*nimages;	//dimension of function value (number of equations)
	double	eps[6];			// precisions for stop-criteria (see manual for more details
	double	*x;				// solution vector. contains values x for f(x)
	MKL_INT	iter1 = 1000;	//maximum number of iterations
	MKL_INT	iter2 = 100;	//maximum number of iterations of calculation of trial-step
	double	rs = 0.0;		//initial step bound
	MKL_INT RCI_Request;	//reverse communication interface parameter
	MKL_INT successful;		//controls of rci cycle
	double *fvec;			//function (f(x)) value vector
	double *fobs;			//The observed value for each function
	double *fjac;			//jacobi matrix
	double	*LW, *UP;		//lower and upper bounds
	MKL_INT iter;			//number of iterations
	MKL_INT st_cr;			//number of stop-criterion
	double r1, r2;			//initial and final residuals
	_TRNSPBC_HANDLE_t handle;	//TR solver handle
	MKL_INT row, col, i, j, k;		//counters
	double h_best, A_best, B_best;	//The best multi-start solution for h, A, and B
	int counter;


	//memory allocation
	x = new double[n];
	fvec = new double[m];
	fobs = new double[m];
	fjac = new double[m*n];
	LW = new double[n];
	UP = new double[n];

	//set precisions for stop-criteria
	for (i = 0; i < 6; i++)
	{
		eps [i] = 0.0000001;
	}

	//set bounds
	for (i = 0; i < npix; i++) {
		LW [i] = h_min;
		UP [i] = h_max;
	}
	for (i=0; i< nimages; i++) {
		LW [npix+i] = 1;
		UP [npix+i]	= 3;
	}
	LW [npix+nimages] = A_min;
	UP [npix+nimages] = A_max;
	LW [npix+nimages+1] = B_min;
	UP [npix+nimages+1] = B_max;

	//Set the initial guess
	for (i = 0; i < npix; i++) x[i] = H[i];
	for (i = 0; i < nimages; i++) x[npix+i] = 1;
	x[npix+nimages] = A;
	x[npix+nimages+1] = B;

	//set initial values
	for (i = 0; i < m; i++)
		fvec [i] = 0.0;
	for (i = 0; i < m*n; i++)
		fjac [i] = 0.0;

	//Set fobs
	counter = 0;
	for (i = 0; i < npix; i++) {
		row = getRow(pixels[i], ncols);
		col = getCol(pixels[i], ncols);
		for (j = 0; j < nimages; j++) {
			fobs[counter] = Images[getIndex(row, col, j, nrows, ncols)];
			counter++;
		}
	}



	//***********************Solve for the best height, A, and B********************************

	/* initialize solver (allocate memory, set initial values)
		handle	in/out:	TR solver handle
		n       in:     number of function variables
		m       in:     dimension of function value
		x       in:     solution vector. contains values x for f(x)
		LW		in:		lower bound
		UP		in:		upper bound
		eps     in:     precisions for stop-criteria
		iter1   in:     maximum number of iterations
		iter2   in:     maximum number of iterations of calculation of trial-step
		rs      in:     initial step bound */
	if (dtrnlspbc_init (&handle, &n, &m, x, LW, UP, eps, &iter1, &iter2, &rs) != TR_SUCCESS)
	{
		// Exit if unsuccessful
		cout << "| error in dtrnlspbc_init\n";
		MKL_FreeBuffers();
		return 1;
	}

	//set initial rci cycle variables
	RCI_Request = 0;
	successful = 0;

	//rci cycle
	while (successful == 0)
	{
		/* call tr solver
			handle		in/out:	tr solver handle
			fvec		in:     vector
			fjac		in:     jacobi matrix
			RCI_request in/out:	return number which denote next step for performing */
		if (dtrnlspbc_solve (&handle, fvec, fjac, &RCI_Request) != TR_SUCCESS)
		{
			// Exit if unsuccessful
			cout << "| error in dtrnlspbc_solve\n";
			MKL_FreeBuffers();
			return 1;
		}

		// according with rci_request value we do next step
		if (RCI_Request == -1 ||
			RCI_Request == -2 ||
			RCI_Request == -3 ||
			RCI_Request == -4 ||
			RCI_Request == -5 ||
			RCI_Request == -6)
			// exit rci cycle
			successful = 1;

		if (RCI_Request == 1)
		{
			// recalculate function value
			for (i = 0; i < npix; i++) {
				for (j = 0; j < nimages; j++) {
					fvec[i*nimages + j] = x[npix+nimages]*x[npix+j]*interpolateFLIC(ZGAP, &FLIC[j*FLIC_points], z_gap_spacing, FLIC_points, x[i]) + x[npix+nimages+1] - fobs[i*nimages + j];
				}
			}
		}
		if (RCI_Request == 2)
		{
			// compute jacobi matrix
			//fjac = df1/dx1, df2/dx1, df3/dx1, df4/dx1,..., df1/dx2, df2/dx2, df2/dx2, df2/dx2,...
			for (i = 0; i < npix; i++) {
				for (j = 0; j < nimages; j++) {
					fjac[i*m + i*nimages + j] = x[npix+nimages]*x[npix+j]*interpolateFLIC(ZGAP, &DERIV[j*FLIC_points], z_gap_spacing, FLIC_points, x[i]);
					fjac[(npix+j)*m + i*nimages + j] = x[npix+nimages]*interpolateFLIC(ZGAP, &FLIC[j*FLIC_points], z_gap_spacing, FLIC_points, x[i]);
					fjac[(npix+nimages)*m + i*nimages + j] = x[npix+j]*interpolateFLIC(ZGAP, &FLIC[j*FLIC_points], z_gap_spacing, FLIC_points, x[i]);
					fjac[(npix+nimages+1)*m + i*nimages + j] = 1;
					//fjac[(npix+2)*m + i*nimages + j] = ...;  if dye dipole angle is also fit
				}
			}
		}
	}

	/* get solution statuses
		handle            in:	TR solver handle
		iter              out:	number of iterations
		st_cr             out:	number of stop criterion
		r1                out:	initial residuals
		r2                out:	final residuals */
	if (dtrnlspbc_get (&handle, &iter, &st_cr, &r1, &r2) != TR_SUCCESS)
	{
		//exit if unsuccessful
		cout << "| error in dtrnlspbc_get\n";
		MKL_FreeBuffers();
		return 1;
	}

	// free handle memory
	if (dtrnlspbc_delete (&handle) != TR_SUCCESS)
	{
		//exit if unsuccessful
		cout << "| error in dtrnlspbc_delete\n";
		MKL_FreeBuffers();
		return 1;
	}

	// Release internal MKL memory
	MKL_FreeBuffers();

	//******************End of Solve for the best height, A, and B****************************

	//Copy the best-fit variables to the output variables
	for (i = 0; i < npix; i++) {
		H[i] = x[i];
	}
	for (i = 0; i< nimages;i++) {
		Ai[i] = x[npix+i];
	}
	A = x[npix+nimages];
	B = x[npix+nimages+1];


	/* free allocated memory */
	delete [] x;
	delete [] fvec;
	delete [] fjac;
	delete [] LW;
	delete [] UP;

	return 1;

}

