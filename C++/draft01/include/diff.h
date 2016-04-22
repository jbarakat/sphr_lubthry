/* DIFFERENCE OPERATORS
 *  Generate spatial difference operators for thin-film evolution
 *  on a uniform, radial grid r, where r[m] = r0 + m*h.
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
 *  df  [input]   derivative by finite differences (size = m)
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

/* IMPLEMENTATIONS */

// first derivative
void d1(int M, double h, double* f, double *df){
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

// second derivative
void d2(int M, double h, double* f, double *df){
	int m;
	double hh = h*h;
	for (m = 0; m < M; m++){
		if (m == 0){					// forward differences
			df[m] = (2.0*f[m] - 5.0*f[m+1] + 4.0*f[m+2] - f[m+3])/( hh);
		}
		else if (m == M-1){		// backward differences
			df[m] = (2.0*f[m] - 5.0*f[m-1] + 4.0*f[m-2] - f[m-3])/(-hh);
		}
		else {								// central differences
			df[m] = (f[m+1] - 2.0*f[m] + f[m-1])/(hh);
		}
	}
}

// third derivative
void d3(int M, double h, double* f, double *df){
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

// fourth derivative
void d4(int M, double h, double* f, double *df){
	int m;
	double hhhh = h*h*h*h;
	for (m = 0; m < M; m++){
		if (m == 0 || m == 1){						// forward differences
			df[m] = (3.0*f[m] - 14.0*f[m+1] + 26.0*f[m+2] - 24.0*f[m+3] + 11.0*f[m+4] - 2.0*f[m+5])/( hhhh);
		}
		else if (m == M-1 || m == M-2){		// backward differences
			df[m] = (3.0*f[m] - 14.0*f[m-1] + 26.0*f[m-2] - 24.0*f[m-3] + 11.0*f[m-4] - 2.0*f[m-5])/(-hhhh);
		}
		else {														// central differences
			df[m] = (f[m+2] - 4.0*f[m+1] + 6.0*f[m] - 4.0*f[m-1] + f[m-2])/(hhhh);
		}
	}
}

// modified bessel operator
void b2(int M, double h, double r0, double* f, double *bf){
	int m;
	double r;
	double d1f[M], d2f[M];
	
	// calculate derivatives
	d1(M, h, f, d1f);
	d2(M, h, f, d2f);
	for (m = 0; m < M; m++){
		r     = r0 + m*h;
		bf[m] = -f[m] + d1f[m]/r + d2f[m];
	}
}

// gradient of modified Bessel operator
void b3(int M, double h, double r0, double* f, double *bf){
	int m;
	double r, r2;
	double d1f[M], d2f[M], d3f[M];
	
	// calculate derivatives
	d1(M, h, f, d1f);
	d2(M, h, f, d2f);
	d3(M, h, f, d3f);
	for (m = 0; m < M; m++){
		r     = r0 + m*h;
		r2    = r*r;
		bf[m] = -(1.0 + 1.0/r2)*d1f[m] + d2f[m]/r + d3f[m];
	}
}

// Laplacian of modified Bessel operator
void b4(int M, double h, double r0, double* f, double *bf){
	int m;
	double r, r2, r3;
	double d1f[M], d2f[M], d3f[M], d4f[M];
	
	// calculate derivatives
	d1(M, h, f, d1f);
	d2(M, h, f, d2f);
	d3(M, h, f, d3f);
	d4(M, h, f, d4f);
	for (m = 0; m < M; m++){
		r     = r0 + m*h;
		r2    = r*r;
		bf[m] = (1.0/r3 - 1.0/r)*d1f[m] - (1.0 + 1.0/r2)*d2f[m] + (2.0/r)*d3f[m] + d4f[m];
	}
}


#endif
