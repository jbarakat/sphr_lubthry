/* DIFFERENCE OPERATORS
 *  Generate spatial difference operators for thin-film evolution
 *  on a uniform, radial grid r, where r[m] = m*h.
 *  All difference operators are second-order in the grid spacing
 *  h, i.e. the error is O(h^2).
 *
 * REFERENCES
 *  Moin P, Fundamentals of Engineering Numerical Analysis
 *  
 * PARAMETERS
 *  m   [input]   space index that ranges from 0 to M-1 (radial grid)
 *  M   [input]   number of grid points
 *  h   [input]   grid spacing
 *  f   [input]   function vector (size = M)
 *  df  [input]   derivative by finite differences (size = M)
 */

#ifndef DIFF_H
#define DIFF_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>

/* PROTOTYPES */
void diff_d1(     int, double, double*, double *);
void diff_d1(int, int, double, double*, double &);

void diff_d2(     int, double, double*, double *);
void diff_d2(int, int, double, double*, double &);

void diff_d3(     int, double, double*, double *);
void diff_d3(int, int, double, double*, double &);

void diff_d4(     int, double, double*, double *);
void diff_d4(int, int, double, double*, double &);

void diff_l2(     int, double, double*, double *);
void diff_l2(int, int, double, double*, double &);

void diff_b2(     int, double, double*, double *);

void diff_b3(     int, double, double*, double *);

void diff_b4(     int, double, double*, double *);

/* IMPLEMENTATIONS */

// first derivative
void diff_d1(int M, double h, double* f, double *df){
	int m;
	double h2 = 2.0*h;
	for (m = 0; m < M; m++){
		if (m == 0){					// forward differences
			df[m] = (-3.0*f[m] + 4.0*f[m+1] - f[m+2])/( h2);
		}
		else if (m == M-1){		// backward differences
			df[m] = (-3.0*f[m] + 4.0*f[m-1] - f[m-2])/(-h2);
		}
		else {								// central differences
			df[m] = (f[m+1] - f[m-1])/(h2);
		}
	}
}

void diff_d1(int m, int M, double h, double* f, double &df){
	double h2 = 2.0*h;
	if (m == 0){					// forward differences
		df = (-3.0*f[m] + 4.0*f[m+1] - f[m+2])/( h2);
	}
	else if (m == M-1){		// backward differences
		df = (-3.0*f[m] + 4.0*f[m-1] - f[m-2])/(-h2);
	}
	else {								// central differences
		df = (f[m+1] - f[m-1])/(h2);
	}
}

// second derivative
void diff_d2(int M, double h, double* f, double *df){
	int m;
	double hh = h*h;
	for (m = 0; m < M; m++){
		if (m == 0){					// forward differences
			df[m] = (2.0*f[m] - 5.0*f[m+1] + 4.0*f[m+2] - f[m+3])/(hh);
		}
		else if (m == M-1){		// backward differences
			df[m] = (2.0*f[m] - 5.0*f[m-1] + 4.0*f[m-2] - f[m-3])/(hh);
		}
		else {								// central differences
			df[m] = (f[m+1] - 2.0*f[m] + f[m-1])/(hh);
		}
	}
}

void diff_d2(int m, int M, double h, double* f, double &df){
	double hh = h*h;
	if (m == 0){					// forward differences
		df = (2.0*f[m] - 5.0*f[m+1] + 4.0*f[m+2] - f[m+3])/(hh);
	}
	else if (m == M-1){		// backward differences
		df = (2.0*f[m] - 5.0*f[m-1] + 4.0*f[m-2] - f[m-3])/(hh);
	}
	else {								// central differences
		df = (f[m+1] - 2.0*f[m] + f[m-1])/(hh);
	}
}

// third derivative
void diff_d3(int M, double h, double* f, double *df){
	int m;
	double hhh2 = 2.0*h*h*h;
	for (m = 0; m < M; m++){
		if (m == 0 || m == 1){						// forward differences
			df[m] = (-5.0*f[m] + 18.0*f[m+1] - 24.0*f[m+2] + 14.0*f[m+3] - 3.0*f[m+4])/( hhh2);
		}
		else if (m == M-1 || m == M-2){		// backward differences
			df[m] = (-5.0*f[m] + 18.0*f[m-1] - 24.0*f[m-2] + 14.0*f[m-3] - 3.0*f[m-4])/(-hhh2);
		}
		else {														// central differences
			df[m] = (f[m+2] - 2.0*f[m+1] + 2.0*f[m-1] - f[m-2])/(hhh2);
		}
	}
}

void diff_d3(int m, int M, double h, double* f, double &df){
	double hhh2 = 2.0*h*h*h;
	if (m == 0 || m == 1){						// forward differences
		df = (-5.0*f[m] + 18.0*f[m+1] - 24.0*f[m+2] + 14.0*f[m+3] - 3.0*f[m+4])/( hhh2);
	}
	else if (m == M-1 || m == M-2){		// backward differences
		df = (-5.0*f[m] + 18.0*f[m-1] - 24.0*f[m-2] + 14.0*f[m-3] - 3.0*f[m-4])/(-hhh2);
	}
	else {														// central differences
		df = (f[m+2] - 2.0*f[m+1] + 2.0*f[m-1] - f[m-2])/(hhh2);
	}
}

// fourth derivative
void diff_d4(int M, double h, double* f, double *df){
	int m;
	double hhhh = h*h*h*h;
	for (m = 0; m < M; m++){
		if (m == 0 || m == 1){						// forward differences
			df[m] = (3.0*f[m] - 14.0*f[m+1] + 26.0*f[m+2] - 24.0*f[m+3] + 11.0*f[m+4] - 2.0*f[m+5])/(hhhh);
		}
		else if (m == M-1 || m == M-2){		// backward differences
			df[m] = (3.0*f[m] - 14.0*f[m-1] + 26.0*f[m-2] - 24.0*f[m-3] + 11.0*f[m-4] - 2.0*f[m-5])/(hhhh);
		}
		else {														// central differences
			df[m] = (f[m+2] - 4.0*f[m+1] + 6.0*f[m] - 4.0*f[m-1] + f[m-2])/(hhhh);
		}
	}
}

void diff_d4(int m, int M, double h, double* f, double &df){
	double hhhh = h*h*h*h;
	if (m == 0 || m == 1){						// forward differences
		df = (3.0*f[m] - 14.0*f[m+1] + 26.0*f[m+2] - 24.0*f[m+3] + 11.0*f[m+4] - 2.0*f[m+5])/(hhhh);
	}
	else if (m == M-1 || m == M-2){		// backward differences
		df = (3.0*f[m] - 14.0*f[m-1] + 26.0*f[m-2] - 24.0*f[m-3] + 11.0*f[m-4] - 2.0*f[m-5])/(hhhh);
	}
	else {														// central differences
		df = (f[m+2] - 4.0*f[m+1] + 6.0*f[m] - 4.0*f[m-1] + f[m-2])/(hhhh);
	}
}

// Laplacian
void diff_l2(int M, double h, double* f, double *lf){
	int m;
	double r;
	double d1f[M], d2f[M];
	
	// calculate derivatives
	diff_d1(M, h, f, d1f);
	diff_d2(M, h, f, d2f);

	// r = 0 (apply L'Hopital's rule)
	m = 0;
	lf[m] = 2.0*d2f[m];
	
	// r > 0
	for (m = 1; m < M; m++){
		r     = m*h;
		lf[m] = d1f[m]/r + d2f[m];
	}
}

void diff_l2(int m, int M, double h, double* f, double &lf){
	double r;
	double d1f, d2f;
	
	// calculate derivatives
	diff_d1(m, M, h, f, d1f);
	diff_d2(m, M, h, f, d2f);

	// r = 0 (apply L'Hopital's rule)
	if (m == 0)
		lf = 2.0*d2f;
	else {
		r  = m*h;
		lf = d1f/r + d2f;
	}
}

// modified Bessel operator
void diff_b2(int M, double h, double* f, double *bf){
	int m;
	double r;
	double d1f[M], d2f[M];
	
	// calculate derivatives
	diff_d1(M, h, f, d1f);
	diff_d2(M, h, f, d2f);

	// r = 0 (apply L'Hopital's rule)
	m = 0;
	bf[m] = -f[m] + 2.0*d2f[m];
	
	// r > 0
	for (m = 1; m < M; m++){
		r     = m*h;
		bf[m] = -f[m] + d1f[m]/r + d2f[m];
	}
}

// gradient of modified Bessel operator
void diff_b3(int M, double h, double* f, double *bf){
	int m;
	double r, r2;
	double d1f[M], d2f[M], d3f[M];
	
	// calculate derivatives
	diff_d1(M, h, f, d1f);
	diff_d2(M, h, f, d2f);
	diff_d3(M, h, f, d3f);

	// r = 0 (apply L'Hopital's rule)
	m = 0;
	bf[m] = -d1f[m] + (3.0/2.0)*d3f[m];
	
	// r > 0
	for (m = 1; m < M; m++){
		r     = m*h;
		r2    = r*r;
		bf[m] = -d1f[m] - d1f[m]/r2 + d2f[m]/r + d3f[m];
	}
}

// Laplacian of modified Bessel operator
void diff_b4(int M, double h, double* f, double *bf){
	int m;
	double r, r2, r3;
	double d1f[M], d2f[M], d3f[M], d4f[M];
	
	// calculate derivatives
	diff_d1(M, h, f, d1f);
	diff_d2(M, h, f, d2f);
	diff_d3(M, h, f, d3f);
	diff_d4(M, h, f, d4f);

	// r = 0 (apply L'Hopital's rule)
	m = 0;
	bf[m] = -2.0*d2f[m] + (8.0/3.0)*d4f[m];
	
	// r > 0
	for (m = 1; m < M-1; m++){
		r     = m*h;
		r2    = r*r;
		bf[m] = -d2f[m] + (1.0 - r2)*d1f[m]/r3 - d2f[m]/r2 + 2.0*d3f[m]/r + d4f[m];
	}
}

#endif
