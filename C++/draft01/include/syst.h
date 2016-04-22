/* DIFFERENCE SYSTEM
 *
 * REFERENCES
 *  
 * PARAMETERS
 *  n   [input]   time index that ranges from 0 to N-1
 *  m   [input]   space index that ranges from 0 to M-1 (radial grid)
 *  N   [input]   number of time steps
 *  M   [input]   number of grid points (system size)
 *  k   [input]   time spacing
 *  h   [input]   grid spacing
 *  f   [input]   unknown function (size = M)
 *  g   [input]   auxiliary (known) function (size = M)
 *  p   [input]   parameters (arbitrary size)
 */

#ifndef DIFF_H
#define DIFF_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include "diff.h"

/* PROTOTYPES */

/* IMPLEMENTATIONS */

// auxiliary function g (if necessary)
void syst_g(int M, double h, double k, double *g){
	int n, m;

}

// unknowns
void syst_lhs(int ID, int M, double h, double k, double* p, double *g, double *f, double *Lf){
	int m;
	if (ID == 0){						/* gravitational spreading of thin film
													 * over horizontal, planar substrate */
		double cf  = p[0] ;	// = (gravitational acceleration)/(kinematic viscosity)
		       cf *= 0.5*k;	// trapezoidal rule weight
		double d1f[M], l2f[M];
		double f2, f3;
		
		// difference operations
		diff_d1(M, h, f, d1f);
		diff_l2(M, h, f, l2f);

		// assemble the left-hand side
		for (m = 0; m < M; m++){
			f2 = f[m]*f[m];
			f3 = f2*f[m];
			Lf[m]  = f[m];
			Lf[m] -= cf*f3*l2f[m]/3.0;
			Lf[m] -= cf*f2*d2f[m]*d2f[m];
		}
	}
}

// knowns
void syst_rhs(int ID, int M, double h, double k, double* p, double *g, double *f, double *Rf){
	int m;
	if (ID == 0){						/* gravitational spreading of thin film
													 * over horizontal, planar substrate */
		double cf  = p[0] ;	// = (gravitational acceleration)/(kinematic viscosity)
		       cf *= 0.5*k;	// trapezoidal rule weight
		double d1f[M], l2f[M];
		double f2, f3;
		
		// difference operations
		diff_d1(M, h, f, d1f);
		diff_l2(M, h, f, l2f);

		// assemble the left-hand side
		for (m = 0; m < M; m++){
			f2 = f[m]*f[m];
			f3 = f2*f[m];
			Rf[m]  = f[m];
			Rf[m] += cf*f3*l2f[m]/3.0;
			Rf[m] += cf*f2*d2f[m]*d2f[m];
		}
	}
}


#endif
