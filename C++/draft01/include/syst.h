/* DIFFERENCE SYSTEM
 *
 * REFERENCES
 *  Lopez et al, J Coll Int Sci (1976) - gravitational spreading
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

#ifndef SYST_H
#define SYST_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include "./diff.h"

/* PROTOTYPES */
void syst_eqn(int,      int, double, double, double *, double *, double *, double *);
void syst_lhs(int,      int, double, double, double *, double *, double *);
void syst_lhs(int, int, int, double, double, double *, double *, double &);
void syst_rhs(int,      int, double, double, double *, double *, double *);
void syst_rhs(int, int, int, double, double, double *, double *, double &);

/* IMPLEMENTATIONS */
// Assemble equation F(f1;f0,p) = 0 and calculate Jacobian DF = dF/df1.
//   NOTE: f0 = knowns, f1 = unknowns, g = auxiliary function, p = parameters
void syst_eqn(int ID, int M, double h, double k, double* p, double *f0, double *f1, double* F, double *DF){
	int m, mi, mj;
	double Lf1[M], Rf0[M];
	double df = 0.0001;
	double f1p[M], Lf1p;

	// assemble left- and right-hand sides
	syst_lhs(ID, M, h, k, p, f1, Lf1);
	syst_rhs(ID, M, h, k, p, f0, Rf0);

	// assemble function vector F
	for (m = 0; m < M; m++){
		F[m] = Lf1[m] - Rf0[m];
	}

	// compute Jacobian matrix DF, exploiting sparsity when possible
	if (ID == 0){						/* gravitational spreading of thin film
                					 * over horizontal, planar substrate */
		// LHS highest derivative = 2
		for (mi = 0; mi < M; mi++){			// element of F
			for (mj = 0; mj < M; mj++){		// element of f
				/*------ CALCULATE PARTIAL DERIVATIVES ------*/
				/*------  WITHOUT EXPLOITING SPARSITY  ------*/
				//for (m = 0; m < M; m++){
				//	f1p[m] = f1[m];
				//}
				//f1p[mj] += df;

				//// calculate partial derivative (forward difference scheme)
				//syst_lhs(ID, mi, M, h, k, p, f1p, Lf1p);
				//DF[mi*M + mj] = (Lf1p - Lf1[mi])/df;

				/*------ CALCULATE PARTIAL DERIVATIVES ------*/
				/*------      EXPLOITING SPARSITY      ------*/
				if      (mi == 0){
					// LHS depends on f1[m], f1[m+1], f1[m+2], f1[m+3]
					if (mj == mi || mj == mi + 1 || mj == mi + 2 || mj == mi + 3){
						// perturb f1
						for (m = 0; m < M; m++){
							f1p[m] = f1[m];
						}
						f1p[mj] += df;

						// calculate partial derivative (forward difference scheme)
						syst_lhs(ID, mi, M, h, k, p, f1p, Lf1p);
						DF[mi*M + mj] = (Lf1p - Lf1[mi])/df;
					}
					else
						DF[mi*M + mj] = 0;
				}
				else if (mi == M){
					// LHS depends on f1[m], f1[m-1], f1[m-2], f1[m-3]
					if (mj == mi || mj == mi - 1 || mj == mi - 2 || mj == mi - 3){
						// perturb f1
						for (m = 0; m < M; m++){
							f1p[m] = f1[m];
						}
						f1p[mj] += df;

						// calculate partial derivative (forward difference scheme)
						syst_lhs(ID, mi, M, h, k, p, f1p, Lf1p);
						DF[mi*M + mj] = (Lf1p - Lf1[mi])/df;
					}
					else
						DF[mi*M + mj] = 0;
				}
				else             {
					// LHS depends on f1[m+1], f1[m], f1[m-1]
					if (mj == mi || mj == mi + 1 || mj == mi - 1){
						// perturb f1
						for (m = 0; m < M; m++){
							f1p[m] = f1[m];
						}
						f1p[mj] += df;

						// calculate partial derivative (forward difference scheme)
						syst_lhs(ID, mi, M, h, k, p, f1p, Lf1p);
						DF[mi*M + mj] = (Lf1p - Lf1[mi])/df;
					}
					else
						DF[mi*M + mj] = 0;
				}
			}
		}
	}
}

// Assemble left-hand side operation (unknowns).
void syst_lhs(int ID, int M, double h, double k, double* p, double *f, double *Lf){
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
		for (m = 1; m < M-1; m++){
			f2     = f[m]*f[m];
			f3     = f2*f[m];
			Lf[m]  = f[m];
			Lf[m] -= cf*f3*l2f[m]/3.0;
			Lf[m] -= cf*f2*d1f[m]*d1f[m];
		}

		// lower BC (Neumann condition)
		m = 0;
		Lf[m] = (-3.0*f[m] + 4.0*f[m+1] - f[m+2])/(2.0*h);

		// upper BC (Dirichlet condition)
		m = M-1;
		Lf[m] = f[m];
	}
}

void syst_lhs(int ID, int m, int M, double h, double k, double* p, double *f, double &Lf){
	if (ID == 0){						/* gravitational spreading of thin film
													 * over horizontal, planar substrate */
		double cf  = p[0] ;	// = (gravitational acceleration)/(kinematic viscosity)
		       cf *= 0.5*k;	// trapezoidal rule weight
		double d1f, l2f;
		double f2, f3;
		
		// difference operations
		diff_d1(m, M, h, f, d1f);
		diff_l2(m, M, h, f, l2f);

		// assemble the left-hand side
		if (m == 0)					// Neumann condition
			Lf  = (-3.0*f[m] + 4.0*f[m+1] - f[m+2])/(2.0*h);
		else if (m == M-1)	// Dirichlet condition
			Lf  = f[m];
		else {
			f2  = f[m]*f[m];
			f3  = f2*f[m];
			Lf  = f[m];
			Lf -= cf*f3*l2f/3.0;
			Lf -= cf*f2*d1f*d1f;
		}
	}
}

// Assemble right-hand side operation (knowns).
void syst_rhs(int ID, int M, double h, double k, double* p, double *f, double *Rf){
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

		// assemble the right-hand side
		for (m = 1; m < M-1; m++){
			f2     = f[m]*f[m];
			f3     = f2*f[m];
			Rf[m]  = f[m];
			Rf[m] += cf*f3*l2f[m]/3.0;
			Rf[m] += cf*f2*d1f[m]*d1f[m];
		}

		// lower BC (Neumann condition)
		m = 0;
		Rf[m] = 0;

		// upper BC (Dirichlet condition)
		m = M-1;
		Rf[m] = 0;
	}
}

void syst_rhs(int ID, int m, int M, double h, double k, double* p, double *f, double &Rf){
	if (ID == 0){						/* gravitational spreading of thin film
													 * over horizontal, planar substrate */
		double cf  = p[0] ;	// = (gravitational acceleration)/(kinematic viscosity)
		       cf *= 0.5*k;	// trapezoidal rule weight
		double d1f, l2f;
		double f2, f3;
		
		// difference operations
		diff_d1(m, M, h, f, d1f);
		diff_l2(m, M, h, f, l2f);

		// assemble the right-hand side
		if (m == 0)
			Rf  = 0;
		else if (m == M-1)
			Rf  = 0;
		else {
			f2  = f[m]*f[m];
			f3  = f2*f[m];
			Rf  = f[m];
			Rf += cf*f3*l2f/3.0;
			Rf += cf*f2*d1f*d1f;
		}
	}
}

// auxiliary function g (if necessary)
void syst_g(int M, double h, double k, double *g){
	int n, m;

}

#endif
