/* DIFFERENCE SYSTEM
 *  General parabolic equation:
 *          du
 *   F(u) = -- - L(u) - g = 0
 *          dt
 *  where L is an operator in one space dimension and g is a source term.
 *  In general, u = u(t,x;p) where t is time, x is the spatial coordinate,
 *  and p is a vector of parameters.
 *  Time is discretized using the trapezoidal method (Crank-Nicholson scheme).
 *
 * REFERENCES
 *  Lopez et al, J Coll Int Sci (1976) - gravitational spreading
 *  
 * PARAMETERS
 *  n   [input]   time index that ranges from 0 to N-1
 *  m   [input]   space index that ranges from 0 to M-1
 *  N   [input]   number of time steps
 *  M   [input]   number of grid points (system size)
 *  o   [input]   order of the spatial difference operator L
 *  k   [input]   time spacing
 *  h   [input]   grid spacing
 *  u   [input]   unknown function (size = M)
 *  g   [input]   source function (size = M)
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

// assemble F and DF for all xn
void syst_A(){
	
}

// parabolic operator F(u) at (tn,xm)
// time difference using the trapezoidal method (Crank-Nicholson)
void syst_F(int id, int o, int n, int m, int M, double h, double k,
            double *p, double *g, double *u0, double *u1, double &F){
	double lhs, rhs;

	syst_lhs(id, o, n, m, M, h, k, p,    u1, lhs);	
	syst_rhs(id, o, n, m, M, h, k, p, g, u0, rhs);	
	F   = lhs - rhs;
}

// left-hand side (unknowns) at (tn,xm)
void syst_lhs(int id, int o, int n, int m, int M, double h, double k,
              double *p,            double *u1, double &lhs){
	double cf = 0.5*k;	// trapezoidal weight
	double Lu1;					// spatial difference operator
	
	syst_L(id, o, n+1, m, N, M, h, k, p, u1, Lu1);
	lhs = u1[m] - cf*Lu1[m];
}

// right-hand side operator (knowns) at (tn,xm)
void syst_rhs(int id, int o, int n, int m, int M, double h, double k,
              double *p, double *g, double *u0, double &rhs){
	double cf = 0.5*k;	// trapezoidal weight
	double Lu0;					// spatial difference operator
	
	syst_L(id, o, n  , m, N, M, h, k, p, u0, Lu0);
	rhs = u0[m] + cf*Lu0[m] + g[m];
}

// BOUNDARY CONDITIONS!!!!!!
// BOUNDARY CONDITIONS!!!!!!
// BOUNDARY CONDITIONS!!!!!!
// BOUNDARY CONDITIONS!!!!!!

// spatial difference operator L(u) at (tn,xm)
void syst_L(int id, int o,
            int n, int m, int N, int M, 
            double h, double k, 
						double *p, double *u, double &Lu){
	double cf           ; // coefficient
	double d1u          ;	// first -order spatial difference operators
	double d2u, l2u, b2u;	// second-     
	double d3u, l3u, b3u;	// third -     
	double d4u, l4u, b4u;	// fourth-     

	/* Heat equation in one dimension with constant diffusivity */
	if (id == 0){
		// difference operators
		diff_d2(m, M, h, u, d2u);

		// parameters
		cf = p[0];						// diffusion coefficient (constant)

		// heat equation operator
		Lu = d2u;
	}

	/* Axisymmetric spreading of a thin film over a horizontal,
	 * planar substrate due to gravity. */
	if (id == 1){
		// difference operators
		diff_d1(m, M, h, u, d1u);
		diff_l2(m, M, h, u, l2u);

		// parameters
		cf = p[0];						/* spreading coefficient cf = rho*g/mu, where
													 *  rho = fluid density
													 *  mu  = fluid viscosity
													 *  g   = acceleration due to gravity */
		
		// gravity-spreading operator
		Lu = cf*u[m]*u[m]*(u[m]*l2u/3.0 + d1u*d1u);
	}
}

// partial derivative dF_i/du1_j
void syst_DF(int id, int o, int n, int mi, int mj, int M, double h, double k,
             double *p, double *g, double *u1, double *lhs, double &DF){
	int m;
	double u1p[M];
	double lhsp  ;
	double du1 = 0.0001;

	for (m = 0; m < M; m++)
		u1p[m] = u1[m];
	u1p[mj] += du1;

	// calculate partial derivative by forward differences
	syst_lhs(id, o, n, mi, M, h, k, p,    u1p, lhsp);	
	DF = (lhsp - lhs[mi])/du1;
}

// jacobian matrix DF = dF/du1
void syst_DF(bool sprs, int id, int o, int n, int M, double h, double k,
             double *p, double *g, double *u1, double *lhs, double *DF){
	int mi, mj;
	double DFij;
		
	if (sprs == false){
		/*------ CALCULATE PARTIAL DERIVATIVES ------*/
		/*------  WITHOUT EXPLOITING SPARSITY  ------*/
		for (mi = 0; mi < M; mi++){			// element of F
			for (mj = 0; mj < M; mj++){		// element of u
				syst_DF(id, o, n, mi, mj, M, h, k, p, g, u1, lhs, DFij);
				DF[mi*M + mj] = DFij;
			}
		}
	}
	
	if (sprs == true){
		/*------ CALCULATE PARTIAL DERIVATIVES ------*/
		/*------      EXPLOITING SPARSITY      ------*/
		for (mi = 0; mi < M; mi++){
			for (mj = 0; mj < M; mj++){
				if   ((o == 2) &&
					   (  (mi == 0   && mj < mi+4)	// forward differences
					   || (mi == M-1 && mj > mi-4)	// backward differences
					   || (mi > 0    && mi < M-1 		// central differences
						     && mj < mi+2 && mj > mi-2)
						 ))
					|| ((o == 3) &&
					   (  (mi == 0   && mj < mi+5) 	// forward  differences
					   || (mi == M-1 && mj > mi-5)	// backward differences
					   || (mi > 0    && mi < M-1 
						 && mj < mi+2 && mj > mi-2)		// central differences
						 ))
					|| ((o == 4) &&
						 (  (mi == 0   && mj < mi+6)	// forward differences
					   || (mi == 1   && mj < mi+6)
					   || (mi == M-1 && mj > mi-6)	// backward differences
					   || (mi == M-2 && mj > mi-6)
					   || (mi > 1    && mi < M-2 		// central differences
						     && mj < mi+3 && mj > mi-3)
						 ))
						syst_DF(id, o, n, mi, mj, M, h, k, p, g, u1, lhs, DFij)
					else
						DFij = 0.0;
				DF[mi*M + mj] = DFij;
			}
		}
	}
}

// START FROM HERE








/////////////////////////////////// OLD WORK BELOW




// Assemble equation F(f1;f0,p) = 0 and calculate Jacobian DF = dF/df1.
//   NOTE: f0 = knowns, f1 = unknowns, g = auxiliary function, p = parameters
void syst_eqn(int id, int M, double h, double k, double* p, double *f0, double *f1, double* F, double *DF){
	int m, mi, mj;
	double Lf1[M], Rf0[M];
	double df = 0.0001;
	double f1p[M], Lf1p;
	int der = 0;

	// assemble left- and right-hand sides
	syst_lhs(id, M, h, k, p, f1, Lf1);
	syst_rhs(id, M, h, k, p, f0, Rf0);

	// assemble function vector F
	for (m = 0; m < M; m++){
		F[m] = Lf1[m] - Rf0[m];
	}

	// determine highest derivative on LHS operator
	if (id == 0)				// heat equation 
		der = 2;
	else if (id == 1)		// evolution equation for gravitational spreading
		der = 2;

	// compute Jacobian matrix DF, exploiting sparsity when possible
	if (der == 2){ // LHS highest derivative = 2
		for (mi = 0; mi < M; mi++){			// element of F
			for (mj = 0; mj < M; mj++){		// element of f
				/*------ CALCULATE PARTIAL DERIVATIVES ------*/
				/*------  WITHOUT EXPLOITING SPARSITY  ------*/
				//for (m = 0; m < M; m++){
				//	f1p[m] = f1[m];
				//}
				//f1p[mj] += df;

				//// calculate partial derivative (forward difference scheme)
				//syst_lhs(id, mi, M, h, k, p, f1p, Lf1p);
				//DF[mi*M + mj] = (Lf1p - Lf1[mi])/df;

				/*------ CALCULATE PARTIAL DERIVATIVES -------*/
				/*------      EXPLOITING SPARSITY      -------*/
				/*- (REQUIRES SPECIFIC KNOWLEDGE OF THE PDE) -*/
				if      (mi == 0){
					// LHS depends on f1[m], f1[m+1], f1[m+2], f1[m+3]
					// unless BC is applied at this grid point
					if (mj == mi || mj == mi + 1 || mj == mi + 2 || mj == mi + 3){
						// perturb f1
						for (m = 0; m < M; m++){
							f1p[m] = f1[m];
						}
						f1p[mj] += df;

						// calculate partial derivative (forward difference scheme)
						syst_lhs(id, mi, M, h, k, p, f1p, Lf1p);
						DF[mi*M + mj] = (Lf1p - Lf1[mi])/df;
					}
					else
						DF[mi*M + mj] = 0;
				}
				else if (mi == M){
					// LHS depends on f1[m], f1[m-1], f1[m-2], f1[m-3]
					// unless BC is applied at this grid point
					if (mj == mi || mj == mi - 1 || mj == mi - 2 || mj == mi - 3){
						// perturb f1
						for (m = 0; m < M; m++){
							f1p[m] = f1[m];
						}
						f1p[mj] += df;

						// calculate partial derivative (forward difference scheme)
						syst_lhs(id, mi, M, h, k, p, f1p, Lf1p);
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
						syst_lhs(id, mi, M, h, k, p, f1p, Lf1p);
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
void syst_lhs(int id, int M, double h, double k, double* p, double *f, double *Lf){
	int m;
	double Lfm;
	
	for (m = 0; m < M; m++){
		syst_lhs(id, m, M, h, k, p, f, Lfm);
		Lf[m] = Lfm;
	}
}

void syst_lhs(int id, int m, int M, double h, double k, double* p, double *f, double &Lf){
	if (id == 1){						/* gravitational spreading of thin film
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
void syst_rhs(int id, int M, double h, double k, double* p, double *f, double *Rf){
	int m;
	double Rfm;
	
	for (m = 0; m < M; m++){
		syst_rhs(id, m, M, h, k, p, f, Rfm);
		Rf[m] = Rfm;
	}
}

void syst_rhs(int id, int m, int M, double h, double k, double* p, double *f, double &Rf){
	if (id == 1){						/* gravitational spreading of thin film
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
		else if (m == M-1){
			Rf  = 0;
			Rf  = 0.0001; // NUMERICAL REASONS, DELETE LATER?
		}
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
