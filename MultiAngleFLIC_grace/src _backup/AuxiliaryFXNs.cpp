/*
 * AuxiliarFXNs.cpp
 *
 *  Created on: 2010-09-09
 *      Author: matt
 */

#include "MultiAngleFLIC.h"

using namespace std;
using namespace Magick;

/*
void generate_image_seq_path(string& path, string image_path, string image_name, string file_ext, int image_num) {
	if (image_num < 10)
		path = image_path + image_name + "00" + IntToStr( image_num ) + file_ext;
	else if (image_num < 100)
		path = image_path + image_name + "0" + IntToStr( image_num ) + file_ext;
	else path = image_path + image_name + IntToStr( image_num ) + file_ext;
}
*/

void generateAveragedRefSeries(string image_path, string base_name, int images_per_set, int nsets) {

	string file_path;
	int nrows, ncols;
	int i, j, k;

	file_path = image_path + base_name + "00_00.tif";
	getImageSize(file_path, nrows, ncols);

	cout << "Image size: " << nrows << " x " << ncols << "\n";

	int* ref_averages = new int[images_per_set*nrows*ncols];
	int* nrefs = new int[images_per_set*nrows*ncols];
	int* image = new int[nrows*ncols];
	int* thresh = new int[nrows*ncols];

	//Initialize ref_averages and nrefs
	for (i = 0; i < images_per_set*nrows*ncols; i++) {
		ref_averages[i] = 0.0;
		nrefs[i] = 0;
	}

	//Load the images
	for (i = 0; i < nsets; i++) {

		cout << "Loading image set: " << i << "\n";

		//Get the file path for the threshold image
		if (i < 10)
			file_path = image_path + base_name + "0" + IntToStr(i) + "_T.tif";
		else
			file_path = image_path + base_name + IntToStr(i) + "_T.tif";

		//Open the threshold image
		if (loadImage(thresh, file_path, nrows, ncols) != 1) {
			cout << "\nAborting Program\n";
			delete [] image;
			delete [] thresh;
			delete [] ref_averages;
			delete [] nrefs;
		}

		//Loop through each image in the set and add to the running reference averages
		for (j = 0; j < images_per_set; j++) {

			cout << "Loading image: " << j << "\n";

			//Get the file path of the image to open
			if (i < 10) {
				if (j < 10)
					file_path = image_path + base_name + "0" + IntToStr(i) + "_0" + IntToStr(j) + ".tif";
				else
					file_path = image_path + base_name + "0" + IntToStr(i) + "_" + IntToStr(j) + ".tif";
			}
			else {
				if (j < 10)
					file_path = image_path + base_name + IntToStr(i) + "_0" + IntToStr(j) + ".tif";
				else
					file_path = image_path + base_name + IntToStr(i) + "_" + IntToStr(j) + ".tif";
			}

			//Open the image
			if (loadImage(image, file_path, nrows, ncols) != 1) {
				cout << "\nAborting Program\n";
				delete [] image;
				delete [] thresh;
				delete [] ref_averages;
				delete [] nrefs;
			}

			//Add to the running reference averages
			for (k = 0; k < nrows*ncols; k++) {
				if (thresh[k] != 0) {  //add it
					ref_averages[j*nrows*ncols + k] = ref_averages[j*nrows*ncols + k] + image[k];
					nrefs[j*nrows*ncols + k]++;
				}
			}
		}
	}

	//Normalize the averages
	for (i = 0; i < images_per_set*nrows*ncols; i++) {
		if (nrefs[i] != 0)
			ref_averages[i] = round(ref_averages[i]/nrefs[i]);
	}

	//Save the images
	for (i = 0; i < images_per_set; i++) {

		//Create the file path
		if (i < 10)
			file_path = image_path + base_name + "_Ave_0" + IntToStr(i) + ".pgm";
		else
			file_path = image_path + base_name + "_Ave_" + IntToStr(i) + ".pgm";

		//Create the image
		writeImagePGM(&ref_averages[i*nrows*ncols], ncols, nrows, file_path);
	}
}

void pixelSeriesGenerator(double* values, int nvals, string image_path, string image_name) {
	int i;
	string path;
	Image new_image( Geometry(1, 1), "white");
	new_image.magick("TIFF");
	new_image.modifyImage();

	for (i = 0; i < nvals; i++) {
		generate_image_seq_path(path, image_path, image_name, ".tif", i);
		new_image.pixelColor(0,0,Color(values[i],values[i],values[i],0));
		new_image.write(path);
	}
}

void generate_image_seq_path(string& path, string image_path, string image_name, string file_ext, int image_num) {
	if (image_num < 10)
		path = image_path + image_name + "0" + IntToStr( image_num ) + file_ext;
	else path = image_path + image_name + IntToStr( image_num ) + file_ext;
	//path = image_path + image_name + IntToStr( image_num ) + file_ext;
}


/*
void generate_image_seq_path(string& path, string image_path, string image_name, string file_ext, int image_num) {
	path = image_path + image_name + IntToStr( image_num ) + file_ext;
}
*/
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

void writeImage(int* image, int width, int height, string output_file) {

	Image new_image( Geometry(width, height), "white");
	new_image.magick("TIFF");
	new_image.modifyImage();

	int i, j, c;

	for (i = 0; i < width; i++) {
		for (j = 0; j < height; j++) {
			c = image[j*width + i];  //assumes row-major ordering in the image array
			if (c < 0) c = 0;
			new_image.pixelColor(i,j,Color(c,c,c,0));
		}
	}
	new_image.write(output_file);
}

void writeImageJPEG(int* image, int width, int height, string output_file) {

	Image new_image( Geometry(width, height), "white");
	new_image.magick("JPEG");
	new_image.modifyImage();

	int i, j, c;

	for (i = 0; i < width; i++) {
		for (j = 0; j < height; j++) {
			c = image[j*width + i];  //assumes row-major ordering in the image array
			if (c < 0) c = 0;
			new_image.pixelColor(i,j,Color(c,c,c,0));
		}
	}
	new_image.write(output_file);
}

void writeImageJPEG(double* image, int width, int height, string output_file) {

	Image new_image( Geometry(width, height), "white");
	new_image.magick("JPEG");
	new_image.modifyImage();

	int i, j, c;

	for (i = 0; i < width; i++) {
		for (j = 0; j < height; j++) {
			c = round(image[j*width + i]);  //assumes row-major ordering in the image array
			if (c < 0) c = 0;
			new_image.pixelColor(i,j,Color(c,c,c,0));
		}
	}
	new_image.write(output_file);
}

void writeImagePGM(int* image, int width, int height, string output_file) {

	Image new_image( Geometry(width, height), "white");
	new_image.magick("PGM");
	new_image.modifyImage();

	int i, j, c;

	for (i = 0; i < width; i++) {
		for (j = 0; j < height; j++) {
			c = image[j*width + i];  //assumes row-major ordering in the image array
			if (c < 0) c = 0;
			new_image.pixelColor(i,j,Color(c,c,c,0));
		}
	}
	new_image.write(output_file);
}

void writeImagePGM(double* image, int width, int height, string output_file) {

	Image new_image( Geometry(width, height), "white");
	new_image.magick("PGM");
	new_image.modifyImage();

	int i, j, c;

	for (i = 0; i < width; i++) {
		for (j = 0; j < height; j++) {
			c = round(image[j*width + i]);  //assumes row-major ordering in the image array
			if (c < 0) c = 0;
			new_image.pixelColor(i,j,Color(c,c,c,0));
		}
	}
	new_image.write(output_file);
}

void writeImage(double* image, double multFactor, int nrows, int ncols, string output_file) {

	Image new_image( Geometry(ncols, nrows), "white");
	new_image.magick("TIFF");
	new_image.modifyImage();

	int i, j, c;

	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			c = round(multFactor*image[j + i*ncols]);  //assumes row-major ordering in the image array
			if (c < 0) c = 0;
			new_image.pixelColor(j,i,Color(c,c,c,0));
		}
	}
	new_image.write(output_file);
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

int readINTCSVtoArray(int* array, string file_path) {

	//Open the CSV file and enter contents into a vector
	vector<int> linear;
	ifstream fin( file_path.c_str() );
	if( fin )
	{
		copy( csv_istream_iterator<int>( fin ), csv_istream_iterator<int>(),
				insert_iterator< vector<int> >( linear, linear.begin() ) );
		fin.close();
	}
	else {
		cout << "\nCSV file with path: " << file_path << " could not be opened\n";
		return 0;
	}

	//Copy the contents of the vector into the double array
	for (int i =0; i < linear.size(); i++) array[i] = (int) linear[i];

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

//Create a list of randomly selected pixels from an image sequence in which each pixel in the list
//has an intensity below cutoff in each of the images in the sequence
void randPixelsC(int* pixels, double* Images, int nrows, int ncols, int nimages, int ndesiredPix, double cutoff) {

	int i, j, k, pix, rowNum, colNum, process;
	int num_NN;				//the number of nearest neighbors for a given pixel
	int NN_row [9];			//the row numbers of the nearest neighbors
	int NN_col [9];			//the column numbers of the nearest neighbors
	list<int> unselected_pix;
	pix = 0;
	srand( time(NULL) );

	//Create a list of the possible pixels that can be selected
	for (i = 0; i < nrows*ncols; i++) {

		rowNum = getRow(i, ncols);
		colNum = getCol(i, ncols);

		//Get the nearest neighbor pixels
		num_NN = getNeighbors(NN_row, NN_col, rowNum, colNum, 0, nrows, ncols);

		//Make sure none of the pixels is above the cutoff threshold
		process = 1;
		for (k =0; k < nimages && process == 1; k++) {
			for (j = 0; j < num_NN; j++) {
				if ( Images[getIndex(NN_row[j], NN_col[j], k, nrows, ncols)] > cutoff) {
					process = 0;
					break;
				}
			}
		}

		if (process == 1) unselected_pix.push_back(i);
	}


	//Select some pixels at random
	for (i = 0; i < ndesiredPix; i++) {

		//select a pixel randomly from the list of unselected pixels
		pix = rand() % unselected_pix.size();  //generates a random # between 0 and .size()-1
		pixels[i] = pix;
		unselected_pix.remove(pix);
	}

}


//Create a list of pixels from an image series. The intensity of the pixel's in the list is below the "cutoff" intensity
//in each of the images of the privided sequence
void allPixelsC(int* pixels, double* Images, int nrows, int ncols, int nimages, int& nPix, double cutoff) {

	int i, j, k, row, col, process;


	//Create a list of the possible pixels that can be selected
	nPix = 0;
	for (i = 0; i < nrows*ncols; i++) {

		row = getRow(i, ncols);
		col = getCol(i, ncols);

		//Make sure none of the pixels is above the cutoff threshold
		process = 1;
		for (k =0; k < nimages && process == 1; k++) {
			if ( Images[getIndex(row, col, k, nrows, ncols)] > cutoff) {
				process = 0;
				break;
			}
		}

		if (process == 1) {
			pixels[nPix] = getIndex(row, col, 1, nrows, ncols);
			nPix++;
		}
	}
}

//Create a list of pixels from a binary threshold image.  The pixels in the list are those that have
//a value of the parameter "val."  Val must be either 0 or 1.  Thresh should be stored in row major format
void threshPixels(int* pixels, int* Thresh, int nrows, int ncols, int& nPix, int val) {

	int i, j, row, col;

	//Create a list of the possible pixels that can be selected
	nPix = 0;
	for (i = 0; i < nrows*ncols; i++) {

		row = getRow(i, ncols);
		col = getCol(i, ncols);

		if (Thresh[i] == val) {
			pixels[nPix] = i;
			nPix++;
		}
	}
}



int getNeighbors(int* m, int* n, int row, int col, int dist, int nrows, int ncols) {

	if (dist == 0) {
		m[0] = row;
		n[0] = col;
		return 1;
	}

	int i, j;
	int count = 0;

	i = row+dist;
	if (i >= 0 && i < nrows) {
		for (j = col-dist; j <= col+dist; j++) {
			if (j >= 0 && j < ncols) {
				m[count] = i;
				n[count] = j;
				count++;
			}
		}
	}

	i = row-dist;
	if (i >= 0 && i < nrows) {
		for (j = col-dist; j <= col+dist; j++) {
			if (j >= 0 && j < ncols) {
				m[count] = i;
				n[count] = j;
				count++;
			}
		}
	}

	j = col+dist;
	if (j >= 0 && j < ncols) {
		for (i = row-dist+1; i <= row+dist-1; i++) {
			if (i >= 0 && i < nrows) {
				m[count] = i;
				n[count] = j;
				count++;
			}
		}
	}

	j = col-dist;
	if (j >= 0 && j < ncols) {
		for (i = row-dist+1; i <= row+dist-1; i++) {
			if (i >= 0 && i < nrows) {
				m[count] = i;
				n[count] = j;
				count++;
			}
		}
	}
	return count;
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

double StrToDouble( const string& s ) {

	double result;
	std::istringstream ss( s );
	ss >> result;
	if (!ss) throw std::invalid_argument( "StrToDouble" );
	return result;

}


double mean(double* x, int nvals) {
	double result = 0.0;
	for (int i= 0; i < nvals; i++) {
		result += x[i];
	}
	return result/nvals;
}

double mean(int* x, int nvals) {
	double result = 0.0;
	for (int i= 0; i < nvals; i++) {
		result += x[i];
	}
	return result/nvals;
}

//Function to remove whitespaces from a string
void removeAllWhiteSpaces(string &str)
{
    string temp;
    for (int i = 0; i < str.length(); i++)
        if (str[i] != ' ') temp += str[i];
    str = temp;
}

void standardFLICFitter (void) {
	double Height  [] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95};
	double Flic [] = {1,0.997057,0.993461,0.989219,0.98434,0.978834,0.972712,0.965987,0.958672,0.950781,0.942329,0.933332,0.923808,0.913774,0.903249,0.892252,0.880803,0.868922,0.856631,0.843952,0.830906,0.817516,0.803805,0.789796,0.775512,0.760977,0.746215,0.731249,0.716102,0.7008,0.685364,0.669819,0.654187,0.638492,0.622756,0.607,0.591248,0.57552,0.559837,0.544219,0.528686,0.513257,0.497951,0.482784,0.467775,0.452939,0.438292,0.423849,0.409623,0.395628,0.381876,0.368378,0.355145,0.342187,0.329512,0.317129,0.305045,0.293265,0.281797,0.270643,0.259809,0.249296,0.239108,0.229246,0.219711,0.210502,0.201619,0.193062,0.184827,0.176913,0.169316,0.162033,0.155059,0.148391,0.142023,0.13595,0.130166,0.124665,0.119441,0.114487,0.109796,0.105363,0.101179,0.0972379,0.0935322,0.090055,0.0867992,0.0837577,0.0809236,0.07829,0.0758504,0.0735984,0.0715276,0.0696321,0.0679062,0.0663445};

	//Open the image
	int nr, nc, ii, jj;
	int * image;
	getImageSize("image.tif", nr, nc);
	image = new int[ nr*nc ];
	double * fits = new double[nr*nc];
	double * x1 = new double[nr*nc];
	double * y1 = new double[nr*nc];

	loadImage(image, "image.tif", nr, nc);
	cout << "rows: " << nr << " cols: " << nc;

	double factor = (255.0/MaxRGB)*(0.52/118.0);
	double max1 = 100;
	int numberfit = 0;
	for (ii = 0; ii < nr*nc; ii++) {
		//cout << "image[" << ii << "]: " << image[ii]*factor << " " << factor << "\n";
		if (image[ii]*factor <= 0.52 && image[ii]*factor >= 0.0715276) {
			jj = 0;
			while (Flic[jj] > image[ii]*factor) {
				//cout << image[ii]*factor << "\n";
				jj++;
			}
			fits[numberfit] = Height[jj];
			x1[numberfit] = getCol(ii, nc);
			y1[numberfit] = getRow(ii, nc);
			numberfit++;
			cout << fits[ii] << " ";
		}
	}

	//write out fit
	arrayWriter(fits, numberfit, "fits.csv");

	cout << "\nrows: " << nr << " cols: " << nc;

	int nx1 = floor(233/0.25) + 1;
	int ny1 = floor(209/0.25) + 1;

	double * xnodes1 = new double[nx1];
	double * ynodes1 = new double[ny1];
	double * zgrid1 = new double[nx1*ny1];

	//Allocate xnodes and ynodes
	for(ii = 0; ii < nx1; ii++)
		xnodes1[ii] = ii*0.25;

	for (ii = 0; ii < ny1; ii++)
		ynodes1[ii] = ii*0.25;

	//Smooth the surface
	Surface (zgrid1, x1, y1, fits, numberfit, xnodes1, ynodes1, nx1, ny1, 5, 't', 'g');

	//Write the best fit heights to file
	arrayWriter(zgrid1, nx1*ny1, "fits_grid2.csv");

	delete [] fits;
	delete [] image;
	delete [] xnodes1;
	delete [] ynodes1;
	delete [] zgrid1;

	cout << "\nnx1: " << nx1 << " ny1: " << ny1 << "\n";
	cout << "\nFinished!";
	return;
}

//Simple csv code
//test out some csv opening
/*
ifstream inFile ("Parameters.dat");
string line;
int linenum = 0;
while (getline (inFile, line))
{
    linenum++;
    cout << "\nLine #" << linenum << ":" << endl;
    istringstream linestream(line);
    string item;
    int itemnum = 0;
    while (getline (linestream, item, ','))
    {
        itemnum++;
        cout << "Item #" << itemnum << ": " << item << endl;
    }
}
*/

/*
int DiO_Lam [] = {475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,
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

int DiI_Lam [] = {530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,
		553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,
		579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,
		606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,
		632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,
		658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,
		684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699};
*/
