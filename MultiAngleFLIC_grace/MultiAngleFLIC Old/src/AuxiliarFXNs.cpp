/*
 * AuxiliarFXNs.cpp
 *
 *  Created on: 2010-09-09
 *      Author: matt
 */

#include "MultiAngleFLIC.h"

using namespace std;
using namespace Magick;

double DiO_Lam [] = {475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,
		496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,
		519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,
		542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,
		565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,
		588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,
		611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,
		634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649};

double DiO_Em [] = {3,3.3,3.6,3.9,4.227634,4.334085,4.851135,5.596293,6.610907,7.863134,9.405489,11.28668,
		13.40097,15.98194,19.05524,22.77534,27.03101,31.90448,37.41,43.53095,49.96317,56.63301,63.31234,69.79446,
		75.69918,81.00748,86.12807,90.85423,94.7511,97.44565,99.2824,100,99.68635,98.55293,96.68528,94.12617,
		90.6974,87.44446,84.21765,80.97185,77.59059,74.24023,71.12748,68.00285,64.98753,62.23357,59.8503,57.75217,
		56.1554,54.78437,53.83153,52.85732,51.91874,51.17025,50.49543,49.99406,49.2432,48.51135,47.77474,47.14744,
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

double DiI_Em [] = {2.843224,2.892825,2.926616,2.983924,3.121656,3.309256,3.48109,3.665323,3.954517,4.373206,
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

double DiI_Lam [] = {530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,
		553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,
		579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,
		606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,
		632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,
		658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,
		684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699};

void getImageSize(string file_path, int& nrows, int& ncols) {

	Image my_image(file_path);

	ncols = my_image.columns();
	nrows = my_image.rows();

}

int loadImage(int* image, string file_path, int nrows, int ncols) {

	Color pix;
	Image my_image(file_path);
	int i, j;


	//Make sure that the image has the expected size
	//Ensure that all images in the sequence have the expected size of (nrows, ncols)
	if ( nrows != my_image.rows() || ncols != my_image.columns() ) {
		cout << "\nError:  All images in a series are not the same size.  Program aborted.\n";
		return 0;
	}

	//Load the pixel values into the image array in row-major format
	for (i = 0; i < ncols; i++) {
		for (j = 0; j < nrows; j++) {
			image[j*ncols + i] = (my_image.pixelColor(i,j)).redQuantum();
		}
	}
	return 1;
}



int get_image_size(int first_image_number, int& nrows, int& ncols, string path, string image_name) {		//Create the proper file path to the image

	uint32 w, h;

	//Create the file_path to the first image in the set
	string file_path = path + image_name + IntToStr( first_image_number ) + ".tif";
	cout << "\nThe path: " << file_path;

    //Open the image
	TIFF* tif;
	tif = TIFFOpen(file_path.c_str(), "r");

	//Get the pixel values, assuming that the image successfully opened
	if (tif) {
		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		nrows = h;
		ncols = w;
		TIFFClose(tif);
		return 1;  //Opening an image and getting its size was successful
	}
	else return 0;  //Opening the image failed, so could not get its size
}

//Function to open images and store their values in a row-major array.
//Images must be tiffs. Relies on libtiff to open the tiff images
int importImage(int* Image, string file_path, int nrows, int ncols) {
	int i, j, k, counter;
	uint32 w, h;
	size_t npixels;
	uint32* raster;
	TIFF* tif;

	tif = TIFFOpen(file_path.c_str(), "r");
	cout << "\n" << file_path << " opened\n";

	//Get the pixel values, assuming that the image successfully opened
	if (tif) {
		TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
		TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
		npixels = w * h;

		//Ensure that all images in the sequence have the expected size of (nrows, ncols)
		if (nrows != (int) h && ncols != (int) w) {
			cout << "\nError:  All images in a series are not the same size.  Program aborted.\n";
			TIFFClose(tif);
			return 0;
		}

		//Get the pixel values from the images and store in the temp array
		raster = (uint32*) _TIFFmalloc(npixels * sizeof (uint32));
		if (raster != NULL) {
			if (TIFFReadRGBAImage(tif, w, h, raster, 0)) {
			    //Add the pixel values to the images array in row major format
				counter =0;
			    for (i = nrows-1; i >= 0; i--) {
			    	for (j = 0; j < ncols; j++) {
			    		Image[counter] = (int) TIFFGetR(raster[i*ncols +j]);
			    		counter++;
			    	}
				}
			}
			_TIFFfree(raster);
		}
		TIFFClose(tif);
	}
	else {
		cout << "\nError:  Could not open a specified image.  Program aborted.\n";
		TIFFClose(tif);
		return 0; //Failed to open an image
	}

	return 1;  //Images was successfully opened
}


//Function to open images and store their values in a row-major array.
//Images must be tiffs. Relies on libtiff to open the tiff images
int importImages(int* images, int first_image_number, int nrows, int ncols, int nimages, string path, string image_name) {
	int i, j, k, counter;
	uint32 w, h;
	size_t npixels;
	uint32* raster;
	TIFF* tif;
	counter =0;
	for (k = 0; k < nimages; k++) {
		cout << "\nhere0 ";
		//Create the proper file path to the image
		string file_path = path + image_name + IntToStr( first_image_number + k ) + ".tif";
		cout << "The path: " << file_path;
	    //Open the image
		tif = TIFFOpen(file_path.c_str(), "r");
		cout << "here1 ";
		//Get the pixel values, assuming that the image successfully opened
		if (tif) {
			TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
			TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
			npixels = w * h;

			//Ensure that all images in the sequence have the expected size of (nrows, ncols)
			if (nrows != (int) h && ncols != (int) w) {
				cout << "\nError:  All images in a series are not the same size.  Program aborted.\n";
				TIFFClose(tif);
				return 0;
			}
			cout << "here2 ";
			//Get the pixel values from the images and store in the temp array
			raster = (uint32*) _TIFFmalloc(npixels * sizeof (uint32));
			if (raster != NULL) {
			    if (TIFFReadRGBAImage(tif, w, h, raster, 0)) {
			    	//Add the pixel values to the images array in row major format
			    	for (i = nrows-1; i >= 0; i--) {
			    		for (j = 0; j < ncols; j++) {
			    			images[counter] = (int) TIFFGetR(raster[i*ncols +j]);
			    			counter++;
			    		}
					}
			    }
			    _TIFFfree(raster);
			}
			TIFFClose(tif);
			cout << "here3 ";

		}
		else {
			cout << "\nError:  Could not open a specified image.  Program aborted.\n";
			TIFFClose(tif);
			return 0; //Failed to open an image
		}
		cout << "here 4";
	}
	cout << " here5 ";
	return 1;  //All images were successfully opened
}

//Function for writing an array to a file
int arrayWriter(double* array, int elements, string file_path) {

	ofstream myfile (file_path.c_str());
	if (myfile.is_open()) {
		for (int i =0; i < elements-1; i++) myfile << array[i] << ",";
		myfile << array[elements-1];
		myfile.close();
	}
	else {
		cout << "Unable to open file with path: " << file_path << " for writing";
		return 0;
	}


	return 1;

}

//Function for writing an array to a file
int arrayWriter(int* array, int elements, string file_path) {

	ofstream myfile (file_path.c_str());
	if (myfile.is_open()) {
		for (int i =0; i < elements-1; i++) myfile << array[i] << ",";
		myfile << array[elements-1];
		myfile.close();
	}
	else {
		cout << "Unable to open file with path: " << file_path << " for writing";
		return 0;
	}


	return 1;

}

int readCSVtoArray(double* array, string file_path) {

	//Open the CSV file and enter contents into a vector
	vector<float> linear;
	ifstream fin( file_path.c_str() );
	if( fin )
	{
		copy( csv_istream_iterator<float>( fin ), csv_istream_iterator<float>(),
				insert_iterator< vector<float> >( linear, linear.begin() ) );
		fin.close();
	}
	else {
		cout << "\nCSV file with path: " << file_path << " could not be opened\n";
		return 0;
	}

	//Copy the contents of the vector into the double array
	for (int i =0; i < linear.size(); i++) array[i] = (double) linear[i];

	//return successful
	return 1;
}

//For row-major image ordering
int getIndex(int row, int col, int image, int nrows, int ncols) {
	return col + row*ncols + image*nrows*ncols;
}

//For row-major image ordering
int getCol(int index, int ncols) {
	return index - ncols*floor(index/ncols);
}

//For row-major image ordering
int getRow(int index, int ncols) {
	return floor(index/ncols);
}

//this still needs to be thought out a bit
double interpolateFLIC(double* curveX, double* curveY, double x_spacing, int sizeCurves, double x) {
	if (x <= 0) return curveY[0];
	if (x >= curveX[sizeCurves - 1]) return curveY[sizeCurves - 1];

	int i = floor(x/x_spacing);
	return (x-curveX[i])/x_spacing*(curveY[i+1]-curveY[i]) + curveY[i];
}

double interpolate(double* curveX, double* curveY, double x) {
	double x_spacing = curveX[1] - curveX[0];
	int i = floor(x/x_spacing);
	return (x-curveX[i])/x_spacing*(curveY[i+1]-curveY[i]) + curveY[i];
}

//Create a list of randomly selected pixels from an image
void randPixels(int* pixels, int nrows, int ncols, int npixels) {

	int i, pix;
	list<int> unselected_pix;
	pix = 0;
	srand( time(NULL) );

	//Create a list of the possible pixels that can be selected
	for (i = 0; i  < nrows*ncols; i++)
		unselected_pix.push_back(i);

	//Select some pixels at random
	for (i = 0; i < npixels; i++) {

		//select a pixel randomly from the list of unselected pixels
		pix = rand() % unselected_pix.size();  //generates a random # between 0 and .size()-1
		pixels[i] = pix;
		unselected_pix.remove(pix);

	}
}

//Create a list of randomly selected pixels from an image
void randPixelsC(int* pixels, double* Images, int nrows, int ncols, int npixels, double cutoff) {

	int i, pix;
	list<int> unselected_pix;
	pix = 0;
	srand( time(NULL) );

	//Create a list of the possible pixels that can be selected
	for (i = 0; i  < nrows*ncols; i++) {
		if (Images[i] < cutoff) {
			unselected_pix.push_back(i);
		}
	}

	//Select some pixels at random
	for (i = 0; i < npixels; i++) {

		//select a pixel randomly from the list of unselected pixels
		pix = rand() % unselected_pix.size();  //generates a random # between 0 and .size()-1
		pixels[i] = pix;
		unselected_pix.remove(pix);

	}
}

int nearestNeighbors(int* m, int* n, int i, int j, int nrows, int ncols) {

	if (i == 0 && j == 0) {
		m[0] = i;
		m[1] = i;
		m[2] = i+1;
		m[3] = i+1;
		n[0] = j;
		n[1] = j+1;
		n[2] = j+1;
		n[3] = j;
		return 4;
	}

	if (i == nrows-1 && j == 0) {
		m[0] = i;
		m[1] = i-1;
		m[2] = i-1;
		m[3] = i;
		n[0] = j;
		n[1] = j;
		n[2] = j+1;
		n[3] = j+1;
		return 4;
	}

	if (i == 0 && j == ncols-1) {
		m[0] = i;
		m[1] = i+1;
		m[2] = i+1;
		m[3] = i;
		n[0] = j;
		n[1] = j;
		n[2] = j-1;
		n[3] = j-1;
		return 4;
	}

	if (i == nrows-1 && j == ncols-1) {
		m[0] = i;
		m[1] = i-1;
		m[2] = i;
		m[3] = i-1;
		n[0] = j;
		n[1] = j;
		n[2] = j-1;
		n[3] = j-1;
		return 4;
	}

	if (i == 0) {
		m[0] = i;
		m[1] = i;
		m[2] = i+1;
		m[3] = i+1;
		m[4] = i+1;
		m[5] = i;
		n[0] = j;
		n[1] = j+1;
		n[2] = j+1;
		n[3] = j;
		n[4] = j-1;
		n[5] = j-1;
		return 6;
	}

	if (i == nrows-1) {
		m[0] = i;
		m[1] = i-1;
		m[2] = i-1;
		m[3] = i;
		m[4] = i;
		m[5] = i-1;
		n[0] = j;
		n[1] = j;
		n[2] = j+1;
		n[3] = j+1;
		n[4] = j-1;
		n[5] = j-1;
		return 6;
	}

	if (j == 0) {
		m[0] = i;
		m[1] = i;
		m[2] = i-1;
		m[3] = i-1;
		m[4] = i+1;
		m[5] = i+1;
		n[0] = j;
		n[1] = j+1;
		n[2] = j+1;
		n[3] = j;
		n[4] = j;
		n[5] = j+1;
		return 6;
	}

	if (j == ncols-1) {
		m[0] = i;
		m[1] = i-1;
		m[2] = i-1;
		m[3] = i;
		m[4] = i+1;
		m[5] = i+1;
		n[0] = j;
		n[1] = j;
		n[2] = j-1;
		n[3] = j-1;
		n[4] = j-1;
		n[5] = j;
		return 6;
	}

    m[0] = i;
    m[1] = i;
    m[2] = i+1;
    m[3] = i+1;
    m[4] = i+1;
    m[5] = i;
    m[6] = i-1;
    m[7] = i-1;
    m[8] = i-1;
    n[0] = j;
    n[1] = j+1;
    n[2] = j+1;
    n[3] = j;
    n[4] = j-1;
    n[5] = j-1;
    n[6] = j-1;
    n[7] = j;
    n[8] = j+1;

    return 9;
}

//Converts a integer to a string - A similar function could be constructed to convert a
//double to a string
string IntToStr( int n ) {

	ostringstream result;
	result << n;
	return result.str();

}

int StrToInt( const string& s ) {

	int result;
	std::istringstream ss( s );
	ss >> result;
	if (!ss) throw std::invalid_argument( "StrToInt" );
	return result;

}

// Implement complex arcsin
//*********************************Has not been tested********************************************
complex<double> cmplxasin(complex<double> z){
	complex<double> casin;
	for(int n=0; n<=25; n++){
		long fac2n = 1;
		for (int i = 2; i <= 2*n; i++){
			fac2n *= i;
		}
		long facn = 1;
		for (int j = 2; j <= n; j++){
			facn *= j;
		}
		complex<double> fac2N(fac2n,0), facN(facn,0), czwei(2,0), czweiN(2*n,0),czweiN1(2*n+1,0);
		casin += (fac2N/(pow(czwei,czweiN)*pow(facN,czwei)))*pow(z,czweiN1)/(czweiN1);
	}
	return casin;
}



