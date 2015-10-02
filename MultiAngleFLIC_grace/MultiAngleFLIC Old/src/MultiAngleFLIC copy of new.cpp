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
#define CHUNKSIZE 200				//parallel region chunk size

static int nthreads;				//The desired number of threads that will not be utilized by parallel algorithms
static double n3 = 1.33;				//The refractive index of the ambient media / cytoplasm
static double n2 = 1.44;				//The refractive index of the cell membrane
static double ng = 1.33;				//The refractive index of the cell membrane / SiO2 gap
static double n1 = 1.46;				//The refractive index of SiO2
static double n0 = 4.00;				//The refractive index of Si
static double d2 = 4.0;				//The thickness of the cell membrane in units of nm
static double d1 = 500.0;			//The thickness of the SiO2 layer
static double thetaDye = PI/2;		//The dye excitation angle (usually pi/2 for DiI/DiO in cells)
static string Dye ("DiO");					//The membrane dye used.  Currently can be either "DiO" or "DiI"
static double lambdaEx = 488.0; //561;		//The wavelength of the microscope excitation laser
static double lambdaEmMin = 500.0; //580;		//The minimum wavelength of the microscope emission filter
static double lambdaEmMax = 550.0; //640;		//The maximum wavelength of the microscope emission filter
static double z_gap_max = 1000.0;		//The maximum possible gap (glycocalyx) thickness (used for constructing the FLIC curve)
static double z_gap_spacing = 1.0;		//The desired sampling rate of points on the FLIC curve
static int FLIC_points = 			//The number of sampled points on the FLIC curve
		floor(z_gap_max/z_gap_spacing)+1;
static double derivative_step = 0.025;
static int num_images = 8;  //28		//The number of images for each sample
static double thetaLaser [] 		//The laser inclination angle (in radians) for each image
    = {0,5,10,15,20,25,30,35};


static double thetaObjMax = 67.0*DTR;  //Upper limit is equal to the critical angle from going from
									//membrane (n2) to medium (n3); Also must consider that the light with this
									//angle from the membrane can enter into the objective (dependent
									//on its NA) after refracting in the water

static int construct_FLIC_curves = 0;//A flag that indicates if the FLIC curve should be constructed (1) or read from file (!=1)
static int FLIC_image_start_num = 0;		//The first number of the image to be analyzed in a sequence
static int ref_image_start_num = 0;  //The first number of reference images
static int background_image_start_num = 0; //The first number of the background images
static int num_image_sets = 2;		//The number of sets of flic images per sample point
static double A_guess = MaxRGB*0.5;	//The initial guess of the scaling parameter for the height fitting algorithm
static double B_guess = MaxRGB*0.05;	//The initial guess of the background parameter for the height fitting algorithm
static const double H_min = 25.0;			//The minimum possible height in the fitting algorithm
static const double H_max = 200.0;			//The maximum possible height in the fitting algorithm
static const double A_min = 30.0;			//The minimum possible scaling parameter in the fitting algorithm
static const double A_max = MaxRGB*0.8;	//The maximum possible scaling parameter in the fitting algorithm
static const double B_min = 0.0;				//The minimum possible background constant in the fitting algorithm
static const double B_max = MaxRGB*0.4;	//The maximum possible background constant in the fitting algorithm
static const double I_cutoff = MaxRGB;//The maximum pixel intensity.  Pixels with greater intensity are not processed
static const int MS_step = 25;				//The sampling rate for initial height guesses in the multi-start method in the
									//fitting algorithm
static int nPixtoFitGlobally = 1000;			//The number of randomly selected pixels that will be globally fit to calculate the best
									//fit A and B


static string FLIC_image_name ("control06_c1_35_"); //("muc02_35_");					//The base name of the main image files (should not include the image sequence #)
static string background_image_name ("control_back_06_35_"); //("muc_back_02_35_");		//The base name of the background image files
static string ref_image_name ("reference");			//The base name of the reference images
static string camera_image_name ("background");				//The base name of the background image for the reference images
static string thresh_image_name ("control06_c1_35_thresh");	 			//The base name of the thresholded image for the image data.  Only pixels below the thresh (value of 0) will be processed
static string image_file_ext (".tif");						//The file extension on the image files.
static string exp_name ("cont06_35"); //("muc02_35t");									//The name of the experiment (used to name output files)
static string image_path ("/home/weaver/workspace/MultiAngleFLIC/images/9_29_10_muc/processed/");	//The path to the directory containing the image files
static string FLIC_path ("/home/weaver/workspace/MultiAngleFLIC/curves/");	//The path to the directory containing the FLIC curves
static string output_path ("/home/weaver/workspace/MultiAngleFLIC/output/");	//The path to the directory that the output files will be written to
static string ZGAP_name ("ZGAP.csv");
static string FLIC_name ("FLIC.csv");
static string DERIV_name ("DERIV.csv");

static double DiO_Lam [] = {475.0,476.0,477.0,478.0,479.0,480.0,481.0,482.0,483.0,484.0,485.0,486.0,487.0,488.0,489.0,
		490.0,491.0,492.0,493.0,494.0,495.0,496.0,497.0,498.0,499.0,500.0,501.0,502.0,503.0,504.0,505.0,506.0,
		507.0,508.0,509.0,510.0,511.0,512.0,513.0,514.0,515.0,516.0,517.0,518.0,519.0,520.0,521.0,522.0,523.0,
		524.0,525.0,526.0,527.0,528.0,529.0,530.0,531.0,532.0,533.0,534.0,535.0,536.0,537.0,538.0,539.0,540.0,
		541.0,542.0,543.0,544.0,545.0,546.0,547.0,548.0,549.0,550.0,551.0,552.0,553.0,554.0,555.0,556.0,557.0,
		558.0,559.0,560.0,561.0,562.0,563.0,564.0,565.0,566.0,567.0,568.0,569.0,570.0,571.0,572.0,573.0,574.0,
		575.0,576.0,577.0,578.0,579.0,580.0,581.0,582.0,583.0,584.0,585.0,586.0,587.0,588.0,589.0,590.0,591.0,
		592.0,593.0,594.0,595.0,596.0,597.0,598.0,599.0,600.0,601.0,602.0,603.0,604.0,605.0,606.0,607.0,608.0,
		609.0,610.0,611.0,612.0,613.0,614.0,615.0,616.0,617.0,618.0,619.0,620.0,621.0,622.0,623.0,624.0,625.0,
		626.0,627.0,628.0,629.0,630.0,631.0,632.0,633.0,634.0,635.0,636.0,637.0,638.0,639.0,640.0,641.0,642.0,
		643.0,644.0,645.0,646.0,647.0,648.0,649.0};

static double DiO_Em [] = {3.0,3.3,3.6,3.9,4.227634,4.334085,4.851135,5.596293,6.610907,7.863134,9.405489,11.28668,
		13.40097,15.98194,19.05524,22.77534,27.03101,31.90448,37.41,43.53095,49.96317,56.63301,63.31234,69.79446,
		75.69918,81.00748,86.12807,90.85423,94.7511,97.44565,99.2824,100,99.68635,98.55293,96.68528,94.12617,
		90.6974,87.44446,84.21765,80.97185,77.59059,74.24023,71.12748,68.00285,64.98753,62.23357,59.8503,57.75217,
		56.1554,54.78437,53.83153,52.85732,51.91874,51.17025,50.49543,49.99406,49.2432,48.51135,47.77474,47.14744,FLIC_path + DERIV_name
		46.57954,45.94749,45.41285,44.7713,44.14875,43.3171,42.4902,41.54687,40.69621,39.55566,38.5458,37.30545,
		36.26233,35.10277,34.16419,33.09731,32.10408,30.88274,29.78258,28.72995,27.77474,26.8623,25.93561,25.18712,
		24.45527,23.69752,22.84876,22.14376,21.41975,20.88559,20.2407,19.75288,19.13865,18.5496,17.92634,17.40596,
		16.83474,16.35286,15.8351,15.43804,15.05738,14.70595,14.37401,14.01758,13.6393,13.2931,12.98277,12.64013,
		12.30557,11.99002,11.76809,11.48224,11.1769,10.8661,10.63918,10.40109,10.13591,9.806819,9.541405,9.298563,
		9.076156,8.735892,8.468575,8.179874,7.955328,7.688012,7.484852,7.260307,7.031009,6.801236,6.58168,6.355471,
		6.096947,5.870738,5.644529,5.474397,5.29286,5.13247,4.915528,4.720209,4.546513,4.415112,4.272069,4.150648,
		3.983367,3.860283,3.712249,3.598432,3.48105,3.343472,3.235119,3.102768,3.026494,2.960675,2.928359,2.844244,
		2.75894,2.654865,2.572651,2.513485,2.41131,2.355637,2.263206,2.265962,2.195248,2.147963,2.063205,2.018296,
		1.934656,1.875252,1.810621,1.799857,1.775858,1.777284,1.774979};

static int DiO_elements = 175;

static double DiI_Lam [] = {530.0,531.0,532.0,533.0,534.0,535.0,536.0,537.0,538.0,539.0,540.0,541.0,542.0,543.0,544.0,545.0,
		546.0,547.0,548.0,549.0,550.0,551.0,552.0,553.0,554.0,555.0,556.0,557.0,558.0,559.0,560.0,561.0,562.0,563.0,
		564.0,565.0,566.0,567.0,568.0,569.0,570.0,571.0,572.0,573.0,574.0,575.0,576.0,577.0,578.0,579.0,580.0,581.0,
		582.0,583.0,584.0,585.0,586.0,587.0,588.0,589.0,590.0,591.0,592.0,593.0,594.0,595.0,596.0,597.0,598.0,599.0,
		600.0,601.0,602.0,603.0,604.0,605.0,606.0,607.0,608.0,609.0,610.0,611.0,612.0,613.0,614.0,615.0,616.0,617.0,
		618.0,619.0,620.0,621.0,622.0,623.0,624.0,625.0,626.0,627.0,628.0,629.0,630.0,631.0,632.0,633.0,634.0,635.0,
		636.0,637.0,638.0,639.0,640.0,641.0,642.0,643.0,644.0,645.0,646.0,647.0,648.0,649.0,650.0,651.0,652.0,653.0,
		654.0,655.0,656.0,657.0,658.0,659.0,660.0,661.0,662.0,663.0,664.0,665.0,666.0,667.0,668.0,669.0,670.0,671.0,
		672.0,673.0,674.0,675.0,676.0,677.0,678.0,679.0,680.0,681.0,682.0,683.0,684.0,685.0,686.0,687.0,688.0,689.0,
		690.0,691.0,692.0,693.0,694.0,695.0,696.0,697.0,698.0,699.0};

static double DiI_Em [] = {2.843224,2.892825,2.926616,2.983924,3.121656,3.309256,3.48109,3.665323,3.954517,4.373206,
		4.946856,5.624003,6.555801,7.733393,9.166076,10.89194,12.85385,15.17715,17.79982,20.71966,23.91408,27.4473,
		31.34721,35.85917,40.58547,46.00974,51.54562,57.65722,63.4411,69.23384,74.56155,79.36227,83.78211,87.62179,
		91.43047,94.41985,96.86005,98.69796,99.83171,100,99.53056,98.31709,97.02834,95.09743,92.85651,90.35873,
		87.29849,84.34455,81.45704,78.53853,75.72188,72.43578,69.35341,66.09832,63.44553,60.86359,58.66696,56.4659,
		54.6147,53.09566,51.47033,49.907,48.38353,47.17892,46.06289,45.14614,44.32241,43.72453,43.16652,42.68822,
		42.34942,41.95084,41.45261,41.35518,41.25731,41.21745,40.9163,40.71568,40.45483,40.07307,39.72808,39.34057,
		38.91187,38.36581,37.76749,37.06599,36.2458,35.39105,34.67804,33.73029,32.82064,31.88087,30.89725,30,28.87865,
		27.94287,26.99513,26.16918,25.33171,24.52613,23.66076,22.86802,21.84012,20.97564,20.12312,19.48228,18.72055,
		17.99247,17.36625,17.00354,16.49513,16.0713,15.34765,14.82152,14.24048,13.85429,13.46457,13.11648,12.73738,
		12.29938,11.97874,11.62622,11.51506,11.17405,11.05226,10.72498,10.5403,10.25332,10.08503,9.910983,9.777679,
		9.64969,9.488485,9.238707,9.114703,8.984056,8.846324,8.591231,8.328167,8.146589,7.930469,7.841452,7.779008,
		7.618689,7.348095,7.084145,7.035872,6.820195,6.622675,6.357396,6.22365,5.986714,5.813995,5.633304,5.551816,
		5.40744,5.389725,5.300266,5.062445,4.857839,4.722764,4.692648,4.624003,4.446413,4.340567,4.232507,4.236935,
		4.161515,4.008547,3.734278};

static int DiI_elements = 170;

static gsl_interp_accel *spectra_acc;
static gsl_spline *spectra_spline;

static gsl_interp_accel *FLIC_acc;
static gsl_spline *FLIC_spline;

static gsl_interp_accel *DERIV_acc;
static gsl_spline *DERIV_spline;


int main()
{


	//***********************************************The local variables**************************************************
	int i, j, k, counter, flag, num_rows, num_cols, index, row, col, chunk;
	double intensity, max, min;
	int *Images_temp;
	double *Images;
	int* Thresh;
	double *H, *R, *A, *B;
	double Ho, Ro, Ao, Bo;
	int *F;
	double *ZGAP, *FLIC, *DERIV;
	int *fitPixels;
	int nfitPix;
	string file_path;
	time_t start;
	double walltime;

	//************Make sure the number of threads does not exceed the maximum available on the system*********************

	if (nthreads < 1) nthreads = 1;
	if ( nthreads > omp_get_max_threads() ); nthreads = omp_get_max_threads();;
	omp_set_num_threads(nthreads);
	cout <<  "Using " << omp_get_max_threads() << " threads in parallel algorithms!!!" << "\n";

	//***********************************Construt the FLIC curves or read them from file**********************************


    //Allocate memory for the FLIC curves,
    ZGAP = new double[FLIC_points];
	FLIC = new double[FLIC_points*num_images];
	DERIV = new double[FLIC_points*num_images];


	//Convert the angles entered as degrees into radians;
	for (i=0; i<num_images; i++) thetaLaser[i] = thetaLaser[i]*DTR;

	//Generate the FLIC curves if the flag, construct_FLIC_Curves,  = 1, otherwise read in FLIC curves from file
	if (construct_FLIC_curves == 1) {

		//Initialize the interpolator for the membrane dye spectra
		spectra_acc = gsl_interp_accel_alloc ();

		if (Dye == "DiO") {
			spectra_spline = gsl_spline_alloc (gsl_interp_cspline, DiO_elements);
			gsl_spline_init (spectra_spline, DiO_Lam, DiO_Em, DiO_elements);
		}
		else if (Dye == "DiI") {
			spectra_spline = gsl_spline_alloc (gsl_interp_cspline, DiI_elements);
			gsl_spline_init (spectra_spline, DiI_Lam, DiI_Em, DiI_elements);
		}
		else {
			cout << "Error: Dye must have a value of 'DiO' or 'DiI' - Aborting Program";
			return 1;
		}

		//Generate the FLIC curves
		start = time(NULL);

		generateFLICCurves(ZGAP, FLIC, DERIV, num_images);  //construct the curves

		walltime = time(NULL) - start;
		printf("generateFLICCurves kernel wall clock time = %.2f sec\n", walltime);

		//Free the interpolator memory
		gsl_spline_free(spectra_spline);
		gsl_interp_accel_free(spectra_acc);

		//Write the FLIC curves to file
		mkdir("curves", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //Make the directory if not already created
		file_path = FLIC_path + ZGAP_name;
		arrayWriter(ZGAP, FLIC_points, "curves/ZGAP.csv");
		file_path = FLIC_path + FLIC_name;
		arrayWriter(FLIC, FLIC_points*num_images, "curves/FLIC.csv");
		file_path = FLIC_path + DERIV_name;
		arrayWriter(DERIV, FLIC_points*num_images, "curves/DERIV.csv");

		cout << "FLIC_points: " << FLIC_points << "\n";

		return 0;
	}

	else {		//otherwise read in FLIC curves from file
		file_path = "curves/ZGAP.csv";
		if ( !readCSVtoArray(ZGAP, file_path) ) {
			cout << "\nFLIC curves (ZGAP.csv) could not be opened - Aborting program\n";
			return 1;
		}
		file_path = "curves/FLIC.csv";
		if ( !readCSVtoArray(FLIC, file_path) ) {
			cout << "\nFLIC curves (FLIC.csv) could not be opened - Aborting program\n";
			return 1;
		}
		file_path = "curves/DERIV.csv";
		if ( !readCSVtoArray(DERIV, file_path) ) {
			cout << "\nFLIC curves (DERIV.csv) could not be opened - Aborting program\n";
			return 1;
		}
	}

	//Initialize the interpolators for the FLIC curves
	FLIC_acc = gsl_interp_accel_alloc ();
	DERIV_acc = gsl_interp_accel_alloc ();
	FLIC_spline = gsl_spline_alloc (gsl_interp_cspline, FLIC_points);
	DERIV_spline = gsl_spline_alloc (gsl_interp_cspline, FLIC_points);
	gsl_spline_init (FLIC_spline, ZGAP, FLIC, FLIC_points);
	gsl_spline_init (DERIV_spline, ZGAP, DERIV, FLIC_points);


    //*****************************Open all the tiff images in the sequence and store****************************
	//*****************************their pixel values in row-major format in "images"****************************
	//Open the background images and record the average background for each of the main FLIC image

	generate_image_seq_path(file_path, image_path, background_image_name, image_file_ext, background_image_start_num);
	getImageSize(file_path, num_rows, num_cols);
	Images_temp = new int[ num_rows*num_cols*num_images*num_image_sets ];
	vector<double> bg_I (num_images,0);
	intensity = 0;

	for (j = 0; j < num_image_sets; j++) {
		for (i = 0; i < num_images; i++) {
			index = j*num_images + i;
			generate_image_seq_path(file_path, image_path, background_image_name, image_file_ext, background_image_start_num + index);
			if ( loadImage(&Images_temp[index*num_rows*num_cols], file_path, num_rows, num_cols) != 1) return 1;

			//Record the average intensity in each of the background images
			intensity = 0;
			for (k = 0; k < num_rows*num_cols; k++) intensity = intensity + Images_temp[k + index*num_rows*num_cols];
			intensity = intensity/(num_rows*num_cols);
			bg_I[i] = bg_I[i] + intensity;
		}
	}

	for (i = 0; i < num_images; i++)
		bg_I[i] = bg_I[i]/num_image_sets;

	delete [] Images_temp;
	cout << "\nFLIC Background: ";
	for (i = 0; i < bg_I.size(); i++)
		cout << bg_I[i] << " ";
	cout << "\n";


/*
	//Open the background image for the reference slide and correct the reference intensities
	file_path = image_path + camera_image_name + image_file_ext;
	getImageSize(file_path, num_rows, num_cols);
	Images_temp = new int[ num_rows*num_cols ];
	double camera_I = 0;
	file_path = image_path + camera_image_name + image_file_ext;

	if (loadImage(Images_temp, file_path, num_rows, num_cols) != 1) return 1;  //import the images

	for (k = 0; k < num_rows*num_cols; k++) {
		camera_I = camera_I + Images_temp[k];
	}

	camera_I = camera_I/(num_rows*num_cols);
	cout << "\nCamera: ";
	cout << camera_I << " ";
	cout << "\n";

	delete [] Images_temp;



	//Open the reference images and record the average intensity for each reference image
	generate_image_seq_path(file_path, image_path, ref_image_name, image_file_ext, ref_image_start_num);
	getImageSize(file_path, num_rows, num_cols);
	Images_temp = new int[ num_rows*num_cols*num_images*num_image_sets ];
	vector<double> ref_I (num_images,0);

	for (j = 0; j < num_image_sets; j++) {
		for (i = 0; i < num_images; i++) {

			index = j*num_images + i;
			generate_image_seq_path(file_path, image_path, ref_image_name, image_file_ext, ref_image_start_num + index);
			if ( loadImage(&Images_temp[index*num_rows*num_cols], file_path, num_rows, num_cols) != 1) return 1;  //import the images

			//Record the average intensity of the pixels in the reference image
			intensity = 0;
			for (k = 0; k < num_rows*num_cols; k++)
				intensity = intensity + Images_temp[k + index*num_rows*num_cols];
			intensity = intensity/(num_rows*num_cols);
			ref_I[i] = ref_I[i] + intensity;
		}
	}

	for (i = 0; i < num_images; i++)
		ref_I[i] = ref_I[i]/num_image_sets - camera_I;

	//Normalize the reference intensities
	min = *min_element(ref_I.begin(),ref_I.end());
	for (i = 0; i < ref_I.size(); i++) ref_I[i] = ref_I[i]/min;

	delete [] Images_temp;


	//Delete this in the future
	ref_I[0] = 1.0;
	ref_I[1] = 1.0;
	ref_I[2] = 1.0;
	ref_I[3] = 1.0;
	ref_I[4] = 1.0;
	ref_I[5] = 1.0;
	ref_I[6] = 1.0;
	ref_I[7] = 1.0;
	ref_I[8] = 1.2;
	ref_I[9] = 1.35;
	ref_I[10] = 1.35;


	cout << "\nFLIC Reference: ";
	for (i = 0; i < ref_I.size(); i++)
		cout << ref_I[i] << " ";
	cout << "\n";

*/

	//Open the Threshold image
	file_path = image_path + thresh_image_name + image_file_ext;
	getImageSize(file_path, num_rows, num_cols);
	Thresh = new int[ num_rows*num_cols ];
	file_path = image_path + thresh_image_name + image_file_ext;

	if (loadImage(Thresh, file_path, num_rows, num_cols) != 1) return 1;  //import the images

	int nThresh = 0;
	//Make sure the threshold image is binary
	for (k = 0; k < num_rows*num_cols; k++) {
		if (Thresh[k] != 0 && Thresh[k] != MaxRGB) {
			cout << "Error:  Threshold image must be binary. Exiting";
			return 1;
		}
		if (Thresh[k] == 0) {
			nThresh++;
		}
	}

	cout <<"\nnThresh: " << nThresh << "rows x cols: " << num_rows*num_cols << "\n";

	//Open the main FLIC image files
	generate_image_seq_path(file_path, image_path, FLIC_image_name, image_file_ext, FLIC_image_start_num);
	getImageSize(file_path, num_rows, num_cols);

	Images = new double[ num_rows*num_cols*num_images ];
	Images_temp = new int[ num_rows*num_cols*num_images*num_image_sets ];

	for (i = 0; i <  num_rows*num_cols*num_images; i++) Images[i] = 0;

	for (j = 0; j < num_image_sets; j++) {
		for (i = 0; i < num_images; i++) {
			index = j*num_images + i;
			generate_image_seq_path(file_path, image_path, FLIC_image_name, image_file_ext, FLIC_image_start_num+index);
			if (loadImage(&Images_temp[index*num_rows*num_cols], file_path, num_rows, num_cols) != 1) return 1;  //import the images
		}
	}

	//Load the average of the images in each image set into images
	for (i = 0; i < num_images; i++) {
		for (k = 0; k < num_rows*num_cols; k++) {
			Images[k + i*num_rows*num_cols] = 0.0;
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
	for (i = 0; i < num_images; i++) {
		for (k = 0; k < num_rows*num_cols; k++) {
			Images[k + i*num_rows*num_cols] = Images[k + i*num_rows*num_cols]/num_image_sets;
		}
	}

	//Correct the intensities of the main FLIC images
	for (i = 0; i < num_images; i++) {
		for (j = 0; j < num_rows*num_cols; j++) {
			//Images[i*num_rows*num_cols + j] = (Images[i*num_rows*num_cols + j] - bg_I[i])/ref_I[i];
			Images[i*num_rows*num_cols + j] = Images[i*num_rows*num_cols + j] - bg_I[i];
		}
	}


	cout << "\nImages:\n";
	row = 73;
	col = 78;
	for (i = 0; i < num_images; i++)
		cout << i << ": " << Images[getIndex(row, col, i, num_rows, num_cols)] << " " << Images_temp[getIndex(row, col, i, num_rows, num_cols)] << "\n";
	cout << "\n";


    //**********************Calculate the best-estimate fit membrane height at each pixel*************************

    //Allocate the
	//matrix H for storing the best-estimate of the membrane height at each image pixel
    H = new double[ num_rows*num_cols ];
    R = new double[ num_rows*num_cols ];
    A = new double[ num_rows*num_cols ];
    B = new double[ num_rows*num_cols ];

    F = new int[ num_rows*num_cols ];

    fitPixels = new int[num_rows*num_cols];
    threshPixels(fitPixels, Thresh, num_rows, num_cols, nfitPix, 0);

    int dist = 0;
    int edge_dist = 0;
    double maxChange = 100;
    int startRow = 125;
    int startCol = 75;

    if (startRow > startCol) edge_dist = startRow;
    else edge_dist = startCol;
    if ( (num_rows - startRow - 1) > edge_dist)
    	edge_dist = num_rows - startRow - 1;
    if ( (num_cols - startCol - 1) > edge_dist)
    	edge_dist = num_cols - startCol - 1;

    dist = 0;
    flag = 1;
    int constraint = 1;

    while (dist <= edge_dist  && flag != 0) {
        flag = FLIC_fitter(Images, startRow, startCol, fitPixels, nfitPix, dist, num_rows, num_cols, num_images,
        		ZGAP, FLIC, DERIV, H, A, B, R, F, constraint, maxChange);

        cout << "dist: " << dist << "\n";
        dist++;
    }

    double* X = new double[ nfitPix ];
    double* Y = new double[ nfitPix ];
    double* HS = new double[ nfitPix ];

	for (i = 0; i < nfitPix; i++) {
		X[i] = getCol(fitPixels[i], num_cols);
		Y[i] = getRow(fitPixels[i], num_cols);
		HS[i] = H[fitPixels[i]];
	}

    //Write the best fit heights to file
    file_path = output_path + exp_name + "_H.csv";
	arrayWriter(HS, nfitPix, file_path);

    //Write the best fit heights to file
    file_path = output_path + exp_name + "_X.csv";
	arrayWriter(X, nfitPix, file_path);

    //Write the best fit heights to file
    file_path = output_path + exp_name + "_Y.csv";
	arrayWriter(Y, nfitPix, file_path);

	file_path = output_path + exp_name + "_R.csv";
	arrayWriter(R, nfitPix, file_path);

	file_path = output_path + exp_name + "_C.csv";
	arrayWriter(F, nfitPix, file_path);


	//**********************************************Smooth Using Surface*************************************


	double xinc = 0.5;
	double yinc = 0.5;
	MKL_INT xmax = num_cols - 1;
	MKL_INT ymax = num_rows - 1;
	MKL_INT nx = floor(xmax/xinc) + 1;
	MKL_INT ny = floor(ymax/yinc) + 1;


	double *xnodes = new double[nx];
	double *ynodes = new double[ny];
	double *zgrid = new double[nx*ny];

	//parameters for Surface
	double smoothness = 2.5;
	char interp = 't';
	char regularizer = 'g';

	//Allocate xnodes and ynodes
	for(i = 0; i < nx; i++)
		xnodes[i] = i*xinc;

	for (i = 0; i < ny; i++)
		ynodes[i] = i*yinc;

	//Smooth the surface
	Surface (zgrid, X, Y, HS, nfitPix, xnodes, ynodes, nx, ny, smoothness, interp, regularizer);

    //Write the best fit heights to file
	//file_path = output_path + exp_name + "_HS.csv";
	file_path = "output/" + exp_name + "_HS.csv";
	arrayWriter(zgrid, nx*ny, file_path);


	//****************************************When no constraints are imposed*********************************
/*
	allPixelsC(fitPixels, Images, num_rows, num_cols, num_images, nfitPix, I_cutoff);
	for (k = 0; k < nfitPix; k++) {

		//Start a multi-start approach to find the best height
		flag = nlsq_FLIC_HAB_pixel_MS (Images, fitPixels[k], num_rows, num_cols, num_images, ZGAP,
		    FLIC, DERIV, Ho, Ao, Bo, Ro, H_min, H_max, MS_step);

		if(flag == 0)  {
			cout << "\nAborting Program!\n";
			MKL_FreeBuffers();
			return 0;
		}

		H[fitPixels[k]] = Ho;
		R[fitPixels[k]] = Ro;
		A[fitPixels[k]] = Ao;
		B[fitPixels[k]] = Bo;

		cout << "\npixel: " << fitPixels[k] << " H best: " << Ho << " A best: " << Ao << " B best: " << Bo << " R: " << Ro;
	}


    file_path = output_path + exp_name + "_H3.csv";
	arrayWriter(H, num_cols*num_rows, file_path);

	file_path = output_path + exp_name + "_R3.csv";
	arrayWriter(R, num_cols*num_rows, file_path);

*/

    //*********************************Delete dynamically allocated memory****************************************
    delete [] Images;
    delete [] Thresh;
    delete [] ZGAP;
    delete [] FLIC;
    delete [] DERIV;
    delete [] H;
    delete [] A;
    delete [] B;
    delete [] R;
    delete [] F;
    delete [] X;
    delete [] Y;
    delete [] HS;
    delete [] fitPixels;
	delete [] xnodes;
	delete [] ynodes;
	delete [] zgrid;

    //Free Intel MKL memory
    MKL_FreeBuffers();

    //Free the interpolator memory
	gsl_spline_free(FLIC_spline);
	gsl_interp_accel_free(FLIC_acc);
	gsl_spline_free(DERIV_spline);
	gsl_interp_accel_free(DERIV_acc);

	cout << "\nrows: " << num_rows << " cols: " << num_cols << " number of fitted Pixels: " << nfitPix << "\n";
	cout << "nx Surface: " << nx << " ny Surface: " << ny << "\n";

    cout << "\nfinished!!!\n\n";
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

	double P = sin(theta2)*Aout*U*gsl_spline_eval(spectra_spline, lambda, spectra_acc);


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
	         = gsl_integration_workspace_alloc (500);


	gsl_function F;
	//F.function = &P_Em_Fromherz_Simple;
	F.function = &P_Em_Fromherz;
    F.params = &parameters;

	gsl_integration_qag(&F, theta2Min, theta2Max, 0.0, 1e-6, 500, 3, w, &result, &error);

	gsl_integration_workspace_free (w);


	return result;
}


double P_Em_theta_lambda(double theta2Min, double theta2Max, double lambdaMin, double lambdaMax, double dg) {

	double param [] = {theta2Min, theta2Max, dg};
	double result, error;

	gsl_integration_workspace * w
			= gsl_integration_workspace_alloc (500);

	gsl_function F;
	F.function = &P_Em_theta;
	F.params = &param;

	gsl_integration_qag(&F, lambdaMin, lambdaMax, 0.0, 1e-6, 500, 3, w, &result, &error);

	gsl_integration_workspace_free (w);

	return result;
}


void generateFLICCurves(double* ZGAP, double* FLIC, double* DERIV, int nimages) {

	int i, j, counter;
	double Pex;
	double *Pem;
	Pem = new double[FLIC_points];
	double max_I = -1.0;


	//Allocate ZGAP, the z-heights of the membrane that the FLIC intensities will be calculated
	//for
	for (i = 0; i<FLIC_points; i++) {
		ZGAP[i] = i*z_gap_spacing;
	}

	//Calculate Pem for each height

	cout << "\nCalculating Pem for FLIC...\n";
	#pragma omp parallel for private(j)
	for (j=0; j < FLIC_points; j++) {
		Pem[j] = P_Em_theta_lambda(0, thetaObjMax, lambdaEmMin, lambdaEmMax, ZGAP[j]);
	}

	//Loop through all the input parameters and calculate Pex and Pex*Pem

	cout << "Calculating Pex for FLIC...\n";
	#pragma omp parallel for private(i,j)
	for (i = 0; i < nimages; i++) {
		for (j = 0; j < FLIC_points; j++) {
			Pex = P_Ex_Fromherz_MA(thetaLaser[i], lambdaEx, ZGAP[j]);
			FLIC[i*FLIC_points + j] = Pex*Pem[j];
		}
	}

	//Calculate the derivatives of the FLIC curves
	//Calculate Pem for each height - Deriv

	cout << "Calculating Pem for FLIC - Deriv...\n";
	#pragma omp parallel for private(j)
	for (j = 0; j < FLIC_points; j++) {
		Pem[j] = P_Em_theta_lambda(0, thetaObjMax, lambdaEmMin, lambdaEmMax, ZGAP[j]+derivative_step);
	}

	//Loop through all the input parameters and calculate Pex and Pex*Pem - Deriv

	cout << "Calculating Pex for FLIC - Deriv...\n";
	#pragma omp parallel for private(i,j)
	for (i = 0; i < nimages; i++) {
		for (j = 0; j < FLIC_points; j++) {
			Pex = P_Ex_Fromherz_MA(thetaLaser[i], lambdaEx, ZGAP[j]+derivative_step);
			DERIV[i*FLIC_points + j] = (Pex*Pem[j] - FLIC[i*FLIC_points + j])/derivative_step;
		}
	}

	//Normalize the FLIC curves
	for (i = 0; i < nimages*FLIC_points; i++) {
		if (FLIC[i] > max_I) {
			max_I = FLIC[i];
		}
	}
	for (i = 0; i < nimages*FLIC_points; i++) {
		FLIC[i] = FLIC[i]/max_I;
		DERIV[i] = DERIV[i]/max_I;
	}
}



//fjac = df1/dx1, df2/dx1, df3/dx1, df4/dx1, df1/dx2, df2/dx2, df2/dx3, df2/dx4,...
/* 	nonlinear least square problem with boundary constraints

	Variables
	images	in:	TR solver handle
	heights out:     solution vector.  con	if ( (startRow - dist < 0) && (startRow + dist >= nrows) && (startCol - dist < 0) && (startCol + dist >= ncols) && (dist >= 4)) {
		return 1;  //exit successful - All pixels processed
	}tains height of each pixel
	A       out:     solution for scaling factor A
	B       out:     solution for background parameter B
*/


int FLIC_fitter(double* Images, int startRow, int startCol, int* fitPixels, int nfit, int& dist, int nrows, int ncols, int nimages,
		double* ZGAP, double* FLIC, double* DERIV, double* H, double* A, double* B, double* R, int* F, int constraint, double maxChange) {

	int* m;  // Holds the row number for neighbor pixels
	int* n;	 // Holds the column number for neighbor pixels
	int nneighbors;  //The number of neighbors
	int pixel;
	int i, j, flag, row, col, num;
	double Ho, Ao, Bo, Ro;
	static int count = 0;
	int mNN [8];  //Holds the row numbers of the nearest-neighbors
	int nNN [8];  //Holds the col numbers of the nearest-neighbors
	int num_NN;  //The number of nearest neighbors
	double h_LB;
	double h_UB;



	if (dist == 0) {

		//initialize the variables
		for (i = 0; i < nrows*ncols; i++) F[i] = 1;  //indicates if a pixel has been processed

		//Process pixels in the fitPixels list
		for (i = 0; i < nfit; i++) {
			F[ fitPixels[i] ] = 0;  //Indicates that it should be processed
		}

		//Get the pixel number for the first pixel to process
		pixel = getIndex(startRow, startCol, 0, nrows, ncols);

		//Calculate the height of the start pixel
		flag = nlsq_FLIC_HAB_pixel_MS (Images, pixel, nrows, ncols, nimages, ZGAP, FLIC, DERIV, Ho, Ao, Bo, Ro, H_min, H_max, MS_step);

		if(flag == 1)  {  //successful
			H[pixel] = Ho;
			R[pixel] = Ro;
			A[pixel] = Ao;
			B[pixel] = Bo;
			F[pixel] = 1;
			//cout << "\nH: " << Ho << " A: " << Ao << " B: " << Bo << " R: " << Ro << "\n";
		}
		else { //unsuccessful
			return 0;
		}
	}

	//Allocate memory to hold the neighbor pixel locations
	m = new int[8*dist];	if ( (startRow - dist < 0) && (startRow + dist >= nrows) && (startCol - dist < 0) && (startCol + dist >= ncols) && (dist >= 4)) {
		return 1;  //exit successful - All pixels processed
	}
	n = new int[8*dist];

	//Get the pixels that are distance "dist" from the startRow and startCol
	nneighbors = getNeighbors(m, n, startRow, startCol, dist, nrows, ncols);

	//Sort and select for the neighbors that have not been processed yet
	num = 0;
	for (i= 0; i < nneighbors; i++) {
		pixel = getIndex(m[i], n[i], 0, nrows, ncols);
		if (F[pixel] != 1) {
			m[num] = m[i];
			n[num] = n[i];
			num++;
		}
	}
	nneighbors = num;  //The number of unprocessed neighbors

	//Calculate the constrained heights of the neighbor pixels
	for (i = 0; i < nneighbors; i++) {

		//Calculate the maximum allowed range of heights based on the height of the nearest neighbors
		num_NN = getNeighbors(mNN, nNN, m[i], n[i], 1, nrows, ncols);

		h_LB = H_max;
		h_UB = H_min;

		for (j = 0; j < num_NN; j++) {
			pixel = getIndex(mNN[j], nNN[j], 0, nrows, ncols);

			if (F[pixel] == 1) {

				if (H[pixel] - maxChange < h_LB)
					h_LB = H[pixel] - maxChange;

				if (H[pixel] + maxChange > h_UB)
					h_UB = H[pixel] + maxChange;
			}
		}

		if (h_LB < H_min)
			h_LB = H_min;
		if (h_UB > H_max)
			h_UB = H_max;


		//Calculate the constrained height
		pixel = getIndex(m[i], n[i], 0, nrows, ncols);
		flag = nlsq_FLIC_HAB_pixel_MS (Images, pixel, nrows, ncols, nimages, ZGAP, FLIC, DERIV, Ho, Ao, Bo, Ro, h_LB, h_UB, MS_step);

		if(flag == 1)  {  //successful
			H[pixel] = Ho;
			R[pixel] = Ro;
			A[pixel] = Ao;
			B[pixel] = Bo;
			F[pixel] = 1;
			cout << " H: " << Ho << " A: " << Ao << " B: " << Bo << " R: " << Ro << "\n";
		}

		else if (flag == 0) { //Catstrophic error
			delete [] n;
			delete [] m;
			return 0;
		}

		else if (flag == 2)
			F[pixel] = 2;  //Pixel processed but the nlsq solver did not converge to an appropriate solution

		else {
			cout << "Unknown flag returned - Aborting\n";
			delete [] n;
			delete [] m;
			return 0;
		}

	}
	// delete allocated memory
	delete [] m;
	delete [] n;

	// Compute heights for the next level of pixels
	/*
	dist = dist + 1;
	cout << "dist: " << dist;
	if ( (startRow - dist < 0) && (startRow + dist >= nrows) && (startCol - dist < 0) && (startCol + dist >= ncols) && (dist >= 4)) {
		return 1;  //exit successful - All pixels processed
	}

	flag = FLIC_fitter_Constrained (Images, startRow, startCol, dist, nrows, ncols, nimages, ZGAP, FLIC, DERIV, H, A, B, R, F, maxChange);

	if (flag == 0) {
		return 0;  //Error occurred in processing - return unsuccessful
	}
	*/
	//Return successful
	return 1;
}


int nlsq_FLIC_HAB_pixel_MS (double* Images, int pixel, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double& H, double& A, double& B, double& R, double h_LB, double h_UB, double stepSize)
{
	//***********************************Variables, etc*******************************************

	int i;		//counters
	double h_best, A_best, B_best;	//The best multi-start solution for h, A, and B
	double r_min;			//Holds the lowest residual from the multi-start solutions
	double MS_start;
	int flag;
	int iMax;
	int success = 0;

	iMax = ceil( (h_UB - h_LB)/stepSize );

	//Start a multi-start approach to find the best height
	r_min = numeric_limits<double>::infinity();
	for (i = 0; i < iMax; i++) {

		H = i*stepSize + h_LB;

		flag = nlsq_FLIC_HAB_pixel(Images, pixel, nrows, ncols, nimages, ZGAP,
				FLIC, DERIV, H, A, B, R, h_LB, h_UB);

		//Check to see if the solution is better than previous multi-start solutions
		if (flag == 1) {  //Returned successful
			if (R < r_min) {
				h_best = H;
				A_best = A;
				B_best = B;
				r_min = R;
			}
			if (success != 1) success = 1;
		}
		else if (flag == 0) //Catastrophic error
			return 0;
	}

	H = h_best;
	A = A_best;
	B = B_best;
	R = r_min;

	if (success == 1)
		return 1;  //successful
	else
		return 2;  //Never converged to a solution
}


int nlsq_FLIC_HAB_pixel(double* Images, int pixel, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double& H, double& A, double& B, double& R, double h_LB, double h_UB)
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
	int NN_row [9];  	//Storage of the row numbers of nearest neighbor pixels
	int NN_col [9];		//Storage of column numbers of nearest neighbor pixels
	int num_NN;				//The number of neighbors for a given pixel
	MKL_INT row, col, i, j, k;		//counters


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
	LW [0] = h_LB;
	LW [1] = A_min;
	LW [2] = B_min;
	UP [0] = h_UB;
	UP [1] = A_max;
	UP [2] = B_max;


	//Loop through all the pixels and calculate the best estimate of the height for each pixel
	row = getRow(pixel, ncols);
	col = getCol(pixel, ncols);

	//Get the nearest neighbor pixels
	num_NN = getNeighbors(NN_row, NN_col, row, col, 0, nrows, ncols);

	//Compute fobs
	for (i = 0; i < num_images; i++) {
		fobs[i] = 0.0;
		for (j = 0; j < num_NN; j++) {
			fobs[i] = fobs[i] + Images[getIndex(NN_row[j], NN_col[j], i, nrows, ncols)];
		}
		fobs[i] = fobs[i]/num_NN;
	}

	//Set the initial guess
	x[0] = H;
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
		return 0;
	}

	//set initial rci cycle variables
	RCI_Request = 0;
	successful = 0;

	//rci cycle
	while (successful == 0)
	{
		//cout << "here1";
		/* call tr solver
			handle		in/out:	tr solver handle
			fvec		in:     vector
			fjac		in:     jacobi matrix
			RCI_request in/out:	return number which denote next step for performing */
		if (dtrnlspbc_solve (&handle, fvec, fjac, &RCI_Request) != TR_SUCCESS)
		{
			// Exit if unsuccessful
			cout << "| error in dtrnlspbc_solve\n";
			return 0;
		}

		//Exit if the solution is infeasible

		if ( isnan(x[0]) || isnan(x[1]) || isnan(x[2]) ) {

			// Exit if unsuccessful
			cout << "| Solver returned NAN\n";
			return -1;
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
			for (i = 0; i < m; i++) {
				fvec[i] = x[1]*interpolateFLIC(ZGAP, &FLIC[i*FLIC_points], z_gap_spacing, FLIC_points, x[0]) + x[2] - fobs[i];

			}
		}
		if (RCI_Request == 2)
		{
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
		return 0;
	}

	//******************End of Solve for the best height, A, and B****************************

	//Check to see if the solution is better than previous multi-start solutions
	H = x[0];
	A = x[1];
	B = x[2];
	R = r2;

	/* free allocated memory */
	delete [] x;
	delete [] fvec;
	delete [] fjac;
	delete [] fobs;
	delete [] LW;
	delete [] UP;

	// free handle memory
	if (dtrnlspbc_delete (&handle) != TR_SUCCESS)
	{
		//exit if unsuccessful
		cout << "| error in dtrnlspbc_delete\n";
		return 0;
	}

	return 1;  //successful
}


//	#pragma omp parallel shared(n,m,eps,LW,UP,iter1,iter2,rs) private(x,RCI_Request,successful,fvec,fobs,fjac,iter,st_cr,r1,r2,handle,NN_row,NN_col,num_NN,row,col,i,j,k,h_best,A_best,B_best,r_min,r1_min,MS_start)
//	#pragma omp for



int FLIC_fitter_Brute (double* Images, int* pixels, int npix, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double A, double B, double* R)
{
	//***********************************Variables, etc*******************************************
	int i, j, k, row, col;		//iterators
	int num_NN;				//the number of nearest neighbors for a given pixel
	int NN_row [9];			//the row numbers of the nearest neighbors
	int NN_col [9];			//the column numbers of the nearest neighbors
	double *fobs =
		new double[nimages];	//the observed data
	double r;				//the residual
	double r_min;			//the minimum residual found
	double h_best;			//the height for a pixel with the lowest residual

	//Loop through all the pixels and calculate the best estimate of the height for each pixel
	for (k = 0; k < npix; k++) {

		row = getRow(pixels[k], ncols);
		col = getCol(pixels[k], ncols);

		//Get the nearest neighbor pixels
		num_NN = getNeighbors(NN_row, NN_col, row, col, 0, nrows, ncols);

		//Compute fobs
		for (i = 0; i < num_images; i++) {
			fobs[i] = 0.0;
			for (j = 0; j < num_NN; j++) {
				fobs[i] = fobs[i] + Images[getIndex(NN_row[j], NN_col[j], i, nrows, ncols)];
			}
			fobs[i] = fobs[i]/num_NN;
		}

		//Find the best height for each pixel with brute force
		r_min = numeric_limits<double>::infinity();

		for (j = 0; ZGAP[j] < H_max; j++) {

			//Compute the residual
			r = 0.0;
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

	delete [] fobs;
	return 1;
}


int FLIC_fitter_BruteT (double* Images, int* pixels, int npix, int nrows, int ncols, int nimages, double* ZGAP,
		double* FLIC, double* DERIV, double* H, double A, double B, double* R)
{
	//***********************************Variables, etc*******************************************
	int i, j, k, row, col;		//iterators
	int num_NN;				//the number of nearest neighbors for a given pixel
	int NN_row [9];			//the row numbers of the nearest neighbors
	int NN_col [9];			//the column numbers of the nearest neighbors
	double *fobs;			//the observed data
	double r;				//the residual
	double r_min;			//the minimum residual found
	double h_best;			//the height for a pixel with the lowest residual

	//Loop through all the pixels and calculate the best estimate of the height for each pixel

	#pragma omp parallel private(i,j,k,row,col,num_NN,NN_row,NN_col,fobs,r,r_min,h_best)
	{
		fobs = new double[nimages];

		#pragma omp for
		for (k = 0; k < npix; k++) {

			row = getRow(pixels[k], ncols);
			col = getCol(pixels[k], ncols);

			//Get the nearest neighbor pixels
			num_NN = getNeighbors(NN_row, NN_col, row, col, 0, nrows, ncols);

			//Compute fobs
			for (i = 0; i < num_images; i++) {
				fobs[i] = 0.0;
				for (j = 0; j < num_NN; j++) {
					fobs[i] = fobs[i] + Images[getIndex(NN_row[j], NN_col[j], i, nrows, ncols)];
				}
				fobs[i] = fobs[i]/num_NN;
			}

			//Find the best height for each pixel with brute force
			r_min = numeric_limits<double>::infinity();

			for (j = 0; ZGAP[j] < H_max; j++) {

				//Compute the residual
				r = 0.0;
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

		delete [] fobs;

	} //#pragma omp parallel shared() private()

	return 1;
}


int nlsq_FLIC_HAB_global (double* Images, int* pixels, int npix, int nrows, int ncols, int nimages, double* ZGAP,
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
		LW [i] = H_min;
		UP [i] = H_max;
	}
	LW [npix] 	= A_min;
	LW [npix+1] = B_min;
	UP [npix] 	= A_max;
	UP [npix+1] = B_max;

	//Set the initial guess
	for (i = 0; i < npix; i++) x[i] = H[ pixels[i] ];
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

	//******************End of Solve for the best height, A, and B****************************

	//Copy the best-fit variables to the output variables
	for (i = 0; i < npix; i++) {
		H[ pixels[i] ] = x[i];
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
