#include "MultiAngleFLIC.h"

using namespace std;


int Surface(double* zgrid, double* x,double* y,double* z, MKL_INT npoints, double* xnodes,double* ynodes, MKL_INT nx,
		MKL_INT ny, double smoothness, char interp, char regularizer) {


	MKL_INT i, j, itemp, idum;
	MKL_INT info, request, sort;
	MKL_INT ngrid;
	char trans;
	double xmin, xmax, ymin, ymax, dtemp, ddum, val;
	MKL_INT *indx, *indy;
	MKL_INT ind, L;
	double dx, dy, tx, ty, t1, t2;

	double *Acsr, *Acsr_t, *Bcsr;	//For sparse storage of the values of matrix A and B in CSR format
	int *ia, *ia_t, *ib;			//For sparse storage of the row indexes for matrix A and B in CSR format
	int *ja, *ja_t, *jb;			//For sparse storage of the columns of matrix A and B in CSR format
	double *Acoo;				//For sparse storage of the values of matrix A in Coord format
	int *Arowind, *Arowind_t;	//For sparse storage of the rows of matrix A and A' in Coord format
	int *Acolind, *Acolind_t;	//For sparse storage of the cols of matrix A and A' in Coord format

	MKL_INT Anrows, Ancols;						//The number of rows / cols of matrix A
	MKL_INT Anrows_t, Ancols_t;					//The number of rows / cols of transpose(A)
	MKL_INT Bnrows, Bncols;  					//The number of rows / cols of matrix B
	MKL_INT Annz, Bnnz, Annz_interp, Annz_reg;	//The number of non-zeros in A and B
	MKL_INT Annzmax, Bnnzmax;					//The maximum number of non zeros that can be stored in A, B

	double *Arhs, *Brhs;  //The right hand sides for equations represented in the A and B matrix

	double Anorm_interp, Anorm_reg;
	double *sums;
	double *temp;

	cout << "\nSurface algorithm called...";

	//***************************Check the parameters for acceptability**************************************
	if (smoothness < 0.0) {
		cout << "Error in Surface: The smoothness parameter must be a positive, definite scalar\n";
		return 0;
	}
	if (interp !=  't' && interp != 'b') {
		cout << "Error in Surface: The 'iterp' parameter must be either 't' or 'b'\n";
		return 0;
	}
	if (regularizer != 'l' && regularizer != 'g'  && regularizer != 'b') {
		cout << "Error in Surface:  the 'regularizer' parameter must be either 'l', 'g', or 'b'\n";
		return 0;
	}

	//*******************Calculate the grid distances between the x nodes and y nodes*************************
	dx = xnodes[1] - xnodes[0];
	dy = ynodes[1] - ynodes[0];

	//*****************************The number of points in the grid*******************************************
	ngrid = nx*ny;

	//*****************************Verify that the nodes are distinct*****************************************
	for (i = 0; i < nx-1; i++) {
		if (xnodes[i+1] <= xnodes[i] || xnodes[i+1] - xnodes[i] != dx) {
			cout << "ERROR in Surface: xnodes and ynodes must be monotone increasing\n";
			return 0;
		}
	}
	for (i = 0; i < ny-1; i++) {
		if (ynodes[i+1] <= ynodes[i] || ynodes[i+1] - ynodes[i] != dy) {
			cout << "ERROR in Surface: xnodes and ynodes must be monotone increasing\n";
			return 0;
		}
	}

	//***********************Find the max and min x and y values in the data points****************************
	xmin = x[ cblas_idamin(npoints, x, 1) ];
	xmax = x[ cblas_idamax(npoints, x, 1) ];
	ymin = y[ cblas_idamin(npoints, y, 1) ];
	ymax = y[ cblas_idamax(npoints, y, 1) ];


	//*********************Do the boundaries of xnodes and ynodes encompass the data?**************************
	if (xmin < xnodes[0]) {
		cout << "ERROR in Surface: xmin > xnodes[0]\n";
		return 0;
	}
	if (ymin < ynodes[0]) {
		ynodes[0] = ymin;
		cout << "ERROR in Surface: ymin > ynodes[0]\n";
		return 0;
	}
	if (xmax > xnodes[nx-1]) {
		xnodes[nx-1] = xmax;
		cout << "ERROR in Surface: xmax > xnodes[nx-1]\n";
		return 0;
	}
	if (ymax > ynodes[ny-1]) {
		ynodes[ny-1] = ymax;
		cout << "ERROR in Surface: ymax > ynodes[ny-1]\n";
		return 0;
	}


	//*****************************************Allocate memory***************************************************|

	//Create the storage arrays for sparse matrix A and B
	Annzmax = 4*npoints + 10*ngrid;

	Acsr = (double*)MKL_malloc( Annzmax*sizeof(double),128);  		//For CSR storage of A
	ia = (int*)MKL_malloc( (npoints+2*ngrid+1)*sizeof(int),128);	//For CSR storage of A
	ja = (int*)MKL_malloc( Annzmax*sizeof(int),128);				//For CSR storage of A

	Bnnzmax = 5*Annzmax;
	Bcsr = (double*)MKL_malloc( Bnnzmax*sizeof(double),128);  		//For CSR storage of B
	ib = (int*)MKL_malloc( (ngrid+1)*sizeof(int),128); 				//For CSR storage of B
	jb = (int*)MKL_malloc( Bnnzmax*sizeof(int),128);				//For CSR storage of B

	Acoo = (double*)MKL_malloc( Annzmax*sizeof(double),128);  		//For Coord storage of A
	Acolind = (int*)MKL_malloc( Annzmax*sizeof(int),128);			//For Coord storage of A
	Arowind = (int*)MKL_malloc( Annzmax*sizeof(int),128);			//For Coord storage of A

	Acsr_t = (double*)MKL_malloc( Annzmax*sizeof(double),128);  	//For CSR storage of A'
	ia_t = (int*)MKL_malloc( (ngrid+1)*sizeof(int),128);			//For CSR storage of A'
	ja_t = (int*)MKL_malloc( Annzmax*sizeof(int),128);				//For CSR storage of A'

	Acolind_t = (int*)MKL_malloc( Annzmax*sizeof(int),128);			//For Coord storage of A'
	Arowind_t = (int*)MKL_malloc( Annzmax*sizeof(int),128);			//For Coord storage of A'

	//The RHS vectors
	Arhs = (double*)MKL_malloc( (npoints + 2*ngrid)*sizeof(double),128);
	Brhs = (double*)MKL_malloc( (ngrid)*sizeof(double),128);

	//Some other structures that will be used later
	sums = (double*)MKL_malloc( (ngrid)*sizeof(double),128);
	temp = (double*)MKL_malloc( (ngrid)*sizeof(double),128);
	indx = (int*)MKL_malloc( (npoints)*sizeof(int),128);
	indy = (int*)MKL_malloc( (npoints)*sizeof(int),128);


	//***********************Determine which cell in the array each point lies in***********************************
	for (i = 0; i < npoints; i++) {
		if (x[i] >= xnodes[nx-1])
			indx[i] = nx-2;
		else {
			for (j = 0; j < nx-1; j++) {
				if (x[i] >= xnodes[j] && x[i] < xnodes[j+1]) {
					indx[i] = j;
					break;
				}
			}
		}
	}

	for (i = 0; i < npoints; i++) {
		if (y[i] >= ynodes[ny-1])
			indy[i] = ny-2;
		else {
			for (j = 0; j < ny-1; j++) {
				if (y[i] >= ynodes[j] && y[i] < ynodes[j+1]) {
					indy[i] = j;
					break;
				}
			}
		}
	}


	//*****************************Interpolation equations for each data point****************************************
	Annz = 0;			//the total number of non-zeros in A
	Annz_interp = 0;  	//the number of non-zeros in A belonging to interpolation equations

	if (interp == 't') {
		//Linear interpolation inside each triangle
		for (i = 0; i < npoints; i++) {

			ind = indy[i] + ny*indx[i];

			dtemp = (x[i] - xnodes[indx[i]])/dx;
			if (dtemp < 0.0) dtemp = 0.0;
			if (dtemp > 1.0) dtemp = 1.0;
			tx = dtemp;

			dtemp = (y[i] - ynodes[indy[i]])/dy;
			if (dtemp < 0.0) dtemp = 0.0;
			if (dtemp > 1.0) dtemp = 1.0;
			ty = dtemp;

			if (tx > ty)
				L = ny;
			else
				L = 1;

			//t1 is min(tx, ty); t2 = max(tx,ty)
			if (tx > ty) {
				t1 = ty;
				t2 = tx;
			}
			else {
				t1 = tx;
				t2 = ty;
			}

			//Create the interpolation equations
			val = 1 - t2;
			if (val != 0.0) {
				Arowind[Annz] = i;
				Acolind[Annz] = ind;
				Acoo[Annz] = val;
				Annz++;
				Annz_interp++;
			}

			val = t1;
			if (val != 0.0) {
				Arowind[Annz] = i;
				Acolind[Annz] = ind + ny +1;
				Acoo[Annz] = val;
				Annz++;
				Annz_interp++;
			}

			val = t2 - t1;
			if (val != 0.0) {
				Arowind[Annz] = i;
				Acolind[Annz] = ind + L;
				Acoo[Annz] = val;
				Annz++;
				Annz_interp++;
			}

			Arhs[i] = z[i];
		}
	}
	else {
		//Bilinear interpolation in a cell
		for (i = 0; i < npoints; i++) {

			ind = indy[i] + ny*indx[i];

			dtemp = (x[i] - xnodes[indx[i]])/dx;
			if (dtemp < 0.0) dtemp = 0.0;
			if (dtemp > 1.0) dtemp = 1.0;
			tx = dtemp;

			dtemp = (y[i] - ynodes[indy[i]])/dy;
			if (dtemp < 0.0) dtemp = 0.0;
			if (dtemp > 1.0) dtemp = 1.0;
			ty = dtemp;

			//Create the interpolation equations
			val = (1 - tx)*(1 - ty);
			if (val != 0.0) {
				Arowind[Annz] = i;
				Acolind[Annz] = ind;
				Acoo[Annz] = val;
				Annz++;
				Annz_interp++;
			}

			val = (1 - tx)*ty;
			if (val != 0.0) {
				Arowind[Annz] = i;
				Acolind[Annz] = ind + 1;
				Acoo[Annz] = val;
				Annz++;
				Annz_interp++;
			}

			val = tx*(1 - ty);
			if (val != 0.0) {
				Arowind[Annz] = i;
				Acolind[Annz] = ind + ny;
				Acoo[Annz] = val;
				Annz++;
				Annz_interp++;
			}

			val = tx*ty;
			if (val != 0.0) {
				Arowind[Annz] = i;
				Acolind[Annz] = ind + ny + 1;
				Acoo[Annz] = val;
				Annz++;
				Annz_interp++;
			}

			Arhs[i] = z[i];
		}
	}

	//******************************Regularizations equations for each data point******************************
	Annz_reg = 0;  //The number of non zeros in A corresponding to regularizations eqns

	// Build regularizer equations in matrix A
	if (regularizer == 'l') {
		//Laplacian regularizer

		//Equations for interior grid points
		for (i = 1; i < nx-1; i++) {
			for (j = 1; j < ny-1; j++) {

				ind = j + ny*i;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind;
				Acoo[Annz] = -2/(dy*dy) - 2/(dx*dx);
				Annz++;
				Annz_reg++;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind - 1;
				Acoo[Annz] = 1/(dy*dy);
				Annz++;
				Annz_reg++;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind + 1;
				Acoo[Annz] = 1/(dy*dy);
				Annz++;
				Annz_reg++;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind - ny;
				Acoo[Annz] = 1/(dx*dx);
				Annz++;
				Annz_reg++;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind + ny;
				Acoo[Annz] = 1/(dx*dx);
				Annz++;
				Annz_reg++;

			}
		}

		//border region
		for (i = 0; i <= nx-1; i = i+nx-1) {
			for (j = 1; j < ny-1; j++) {

				ind = j + ny*i;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind - 1;
				Acoo[Annz] = 1/(dy*dy);
				Annz++;
				Annz_reg++;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind;
				Acoo[Annz] = -2/(dy*dy);
				Annz++;
				Annz_reg++;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind + 1;
				Acoo[Annz] = 1/(dy*dy);
				Annz++;
				Annz_reg++;
			}
		}

		//border region
		for (j = 0; j <= ny-1; j = j+ny-1) {
			for (i = 1; i < nx-1; i++) {

				ind = j + ny*i;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind - ny;
				Acoo[Annz] = 1/(dx*dx);
				Annz++;
				Annz_reg++;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind;
				Acoo[Annz] = -2/(dx*dx);
				Annz++;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind + ny;
				Acoo[Annz] = 1/(dx*dx);
				Annz++;
				Annz_reg++;
			}
		}

		Anrows = npoints + ngrid;   //the number of rows in A
		Ancols = ngrid;				//the number of columns in A

		//All the laplacian regularization equations should be set to zero
		for (i = npoints; i < Anrows; i++) {
			Arhs[i] = 0.0;
		}
	}

	else if (regularizer == 'g'){
		//Gradient regularizer
		for (i = 0; i < nx; i++) {
			for (j = 1; j < ny-1; j++) {

				ind = j + ny*i;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind - 1;
				Acoo[Annz] = 1/(dy*dy);
				Annz++;
				Annz_reg++;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind;
				Acoo[Annz] = -2/(dy*dy);
				Annz++;
				Annz_reg++;

				Arowind[Annz] = ind + npoints;
				Acolind[Annz] = ind + 1;
				Acoo[Annz] = 1/(dy*dy);
				Annz++;
				Annz_reg++;

			}
		}

		for (i = 1; i < nx-1; i++) {
			for (j = 0; j < ny; j++) {

				ind = j + ny*i;

				Arowind[Annz] = ind + npoints + ngrid;
				Acolind[Annz] = ind - ny;
				Acoo[Annz] = 1/(dx*dx);
				Annz++;
				Annz_reg++;

				Arowind[Annz] = ind + npoints + ngrid;
				Acolind[Annz] = ind;
				Acoo[Annz] = -2/(dx*dx);
				Annz++;
				Annz_reg++;

				Arowind[Annz] = ind + npoints + ngrid;
				Acolind[Annz] = ind + ny;
				Acoo[Annz] = 1/(dx*dx);
				Annz++;
				Annz_reg++;

			}
		}

		Anrows = npoints + 2*ngrid;  //The number of rows of A
		Ancols = ngrid;					//The number of columns of A

		//All the gradient equations should be set to zero
		for (i = npoints; i < Anrows; i++) {
			Arhs[i] = 0.0;
		}
	}


	//****************************Calculate matrix 1-norms***********************************************

	//Calculate the 1-norm of a sub matrix of A corresponding to the interpolation eqns
	for (i = 0; i < ngrid; i++)
		sums[i] = 0.0;

	for (i = 0; i < Annz_interp; i++)
		sums[Acolind[i]] = sums[Acolind[i]] + abs(Acoo[i]);

	Anorm_interp = sums[ cblas_idamax(ngrid, sums, 1) ];


	//Calculate the 1-norm of the sub matrix of A corresponding to the regularization eqns
	for (i = 0; i < ngrid; i++)
		sums[i] = 0.0;

	for (i = Annz_interp; i < Annz; i++)
		sums[Acolind[i]] = sums[Acolind[i]] + abs(Acoo[i]);

	Anorm_reg = sums[ cblas_idamax(ngrid, sums, 1) ];

	cout << "\nAnnz: " << Annz << " " << " " << Anrows << " "  << Anorm_interp << " " << Anorm_reg;


	//****************************Normalize matrix A****************************************************
	//Normalize the values of the regularization equations
	for (i = Annz_interp; i < Annz; i++)
		Acoo[i] = Acoo[i]*smoothness*Anorm_interp/Anorm_reg;


	//******************************Convert A to the CSR storage format*********************************
	MKL_INT job [6];
	job[0] = 1;  //Convert from Coordinate to CSR format
	job[1] = 1;  //One based indexing for the CSR matrix
	job[2] = 0;  //Zero based indexing for the Coordinate matrix
	job[4] = Annz;  //The number of non-zeros
	job[5] = 0;  //Fill in Acsr, ja, and ia

	mkl_dcsrcoo (job, &Anrows, Acsr, ja, ia, &Annz, Acoo, Arowind, Acolind, &info);



	//************************************Compute the transpose of A**************************************

	//Create the transpose of A in Coordinate storage format
	for (i = 0; i < Annz; i++) {
		Acolind_t[i] = Arowind[i];
		Arowind_t[i] = Acolind[i];
	}

	Anrows_t = ngrid;
	Ancols_t = Anrows;

	//Convert to CSR storage format
	mkl_dcsrcoo (job, &Anrows_t, Acsr_t, ja_t, ia_t, &Annz, Acoo, Arowind_t, Acolind_t, &info);


	//**************************************Compute B = A'*A**********************************************
	sort = 0;		//Elements in A are already sorted
	trans = 'n';  	//Do not multiply by the transposes of the supplied matrices
	request = 0;	//Memory is already allocated for matrix B - perform multiplication

	//Compute B = A'*A
	mkl_dcsrmultcsr (&trans, &request, &sort, &Anrows_t, &Ancols_t, &Ancols, Acsr_t, ja_t, ia_t, Acsr, ja, ia,
		Bcsr, jb, ib, &Bnnzmax, &info);

	//Matrix B is ngrid x ngrid
	Bnrows = ngrid;
	Bncols = ngrid;


	//***********************************Reorder Matrix B********************************************************

	//B is a symmetric matrix.  Remove elements from storage that are not in the upper triangle and stored
	//redundantly.  Remember that elements are stored with 1-based indexing, which is required by the sparse
	//solver.  Elements on the diagnol with value of zero must also be stored for the solver.
	Bnnz = 0;
	MKL_INT diag_added;
	for (i = 1; i <= Bnrows; i++) {
		itemp = Bnnz;
		diag_added = 0;
		for (j = ib[i-1]; j < ib[i]; j++) {
			if(jb[j-1] >= i) {
				if (jb[j-1] == i) diag_added = 1;
				if (diag_added == 0) {
					//Add the zero diaganol element before adding the current element
					Bcsr[Bnnz] = 0.0;
					jb[Bnnz] = i;
					Bnnz++;
					diag_added = 1;
				}
				Bcsr[Bnnz] = Bcsr[j-1];
				jb[Bnnz] = jb[j-1];
				Bnnz++;
			}
		}
		ib[i-1] = itemp+1;
	}
	ib[Bnrows] = Bnnz+1;

	cout << "\nBnnz: " << Bnnz << "\n";


	//***************************Calculate the RHS for B: Brhs = A'*Arhs******************************************

	trans = 'n';
	mkl_dcsrgemv(&trans, &Anrows_t, Acsr_t, ia_t, ja_t, Arhs, Brhs);


	//*****************************************Solve the system of equations****************************************

	//Solve the system using Intel MKL's direct sparse solver - PARDISO interface
	cout << "\nSurface - Direct solver called";

	//Matrix descriptors
	MKL_INT n = ngrid;
	MKL_INT mtype = 2; // Real symmetric matrix - Change to -2 if not pos def.
	MKL_INT nrhs = 1; 	// Number of right hand sides.

	// Internal solver memory pointer.  OK with 32 and 64-bit arch
	void *pt[64];

	//Pardiso control parameters.
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;


	//Setup the pardizo control parameters
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}

	iparm[0] = 1; // No solver default
	iparm[1] = 2; // Fill-in reordering from METIS
	iparm[2] = 1; // Numbers of processors = value of OMP_NUM_THREADS
	iparm[7] = 2; // Max numbers of iterative refinement steps
	iparm[9] = 13; // Perturb the pivot elements with 1E-13
	iparm[10] = 1; // Use nonsymmetric permutation and scaling MPS
	iparm[26] = 1;  //check matrix
	maxfct = 1; // Maximum number of numerical factorizations.
	mnum = 1; // Which factorization to use.
	msglvl = 0; // Print statistical information in file
	error = 0; // Initialize error flag

	//Initialize the internal solver memory pointer. This is only
    //necessary for the FIRST call of the PARDISO solver.
	for (i = 0; i < 64; i++) {
		pt[i] = 0;
	}


	//Reordering and Symbolic Factorization. This step also allocates
	//all memory that is necessary for the factorization.
	phase = 11;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Bnrows, Bcsr, ib, jb, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);

	if (error != 0) {
		cout << "\nSurface - ERROR during symbolic factorization: " << error;
		return 11;
	}
	cout << "\nSurface - Reordering completed ... ";

	// Numerical factorization.
	phase = 22;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Bnrows, Bcsr, ib, jb, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);

	if (error != 0) {
		cout << "\nSurface - ERROR during numerical factorization: " <<  error;
		return 22;
	}
	cout << "\nSurface - Factorization completed ... ";

	//Back substitution and iterative refinement.
	phase = 33;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Bnrows, Bcsr, ib, jb, &idum, &nrhs,
		iparm, &msglvl, Brhs, zgrid, &error);

	if (error != 0) {
		cout << "\nSurface - ERROR during solution: " <<  error;
		return 33;
	}

	cout << "\nSurface - Solve completed ...\n";

	//Termination and release of memory.
	phase = -1; // Release internal memory.
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &Bnrows, &ddum, ib, jb, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);


	//**************************************Reorganize zgrid**********************************************
	//Currently zgrid is stored in column-major format.  Convert to row major format in the temp vector
	int row, col;
	for (i = 0; i < ngrid; i++) {
		col = floor(i/ny);
		row = i - ny*col;
		//col = (nx-1) - col;
		temp[col + row*nx] = zgrid[i];
	}

	//Copy temp into zgrid.  zgrid will then be in row-major format.
	for ( i = 0; i < ngrid; i++) zgrid[i] = temp[i];

	//*******************************************Free Memory**********************************************
	//Free any internal intel MKL memory
	MKL_Free_Buffers();

	//Free other memory
	MKL_free(Acoo);
	MKL_free(Arowind);
	MKL_free(Acolind);
	MKL_free(Arowind_t);
	MKL_free(Acolind_t);
	MKL_free(Acsr);
	MKL_free(ia);
	MKL_free(ja);
	MKL_free(Acsr_t);
	MKL_free(ia_t);
	MKL_free(ja_t);
	MKL_free(Arhs);
	MKL_free(Bcsr);
	MKL_free(ib);
	MKL_free(jb);
	MKL_free(Brhs);
	MKL_free(indx);
	MKL_free(indy);
	MKL_free(sums);
	MKL_free(temp);

	//*********************************Algorithm completed.  Return successful****************************
	return 1;
}
