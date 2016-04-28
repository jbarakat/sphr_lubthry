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
 *  n    [input]   time index that ranges from 0 to N-1
 *  m    [input]   space index that ranges from 0 to M-1
 *  N    [input]   number of time steps
 *  M    [input]   number of grid points (system size)
 *  o    [input]   order of the spatial difference operator L
 *  k    [input]   time spacing
 *  h    [input]   grid spacing
 *  u    [input]   unknown function (size = M)
 *  g    [input]   (known) source function (size = M)
 *  p    [input]   parameters (arbitrary size)
 *  Lid  [input]	 indicator for problem type
 *  sprs [input]	 indicator for whether or not to exploit sparsity
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
void syst_F_DF(bool, int, int, int,           int, double, double, double *, double *, double *, double *, double *, double *, double *);
void syst_F   (      int, int, int, int,      int, double, double, double *, double *, double *, double *, double *, double &          );
void syst_F   (      int, int, int,           int, double, double, double *, double *, double *, double *, double *, double *          );
void syst_lhs (      int, int, int, int,      int, double, double, double *,           double *,           double *, double &          );
void syst_rhs (      int, int, int, int,      int, double, double, double *, double *,           double *,           double &          );
void syst_DF  (      int, int, int, int, int, int, double, double, double *,           double *,           double *, double *, double &);
void syst_DF  (bool, int, int, int,           int, double, double, double *,           double *,           double *, double *, double *);
void syst_L   (char, int, int, int, int,      int, double, double, double *,           double *,                     double &          );

/* IMPLEMENTATIONS */

// assemble F and DF for all xm at fixed tn
void syst_F_DF(bool sprs, int Lid, int o, int n, int M, double h, double k,
	             double *p, double *g0, double *g1, double *u0, double *u1, double *F, double *DF){
	int m;
	double Fm, lhsm, rhsm;
	double lhs[M], rhs[M], DFm[M*M];

	// calculate F
	for (m = 0; m < M; m++){
		syst_lhs(Lid, o, n, m, M, h, k, p, g1, u1, lhsm);	
		syst_rhs(Lid, o, n, m, M, h, k, p, g0, u0, rhsm);	
		Fm   = lhsm - rhsm;

		lhs[m] = lhsm;
		rhs[m] = rhsm;
		F  [m] = Fm  ;
	}
	
	// calculate DF
	syst_DF(sprs, Lid, o, n, M, h, k, p, g1, u1, lhs, DFm);
	for (m = 0; m < M*M; m++){
		DF[m] = DFm[m];
	}
}

// parabolic operator F(u) at (tn,xm)
// time difference using the trapezoidal method (Crank-Nicholson)
void syst_F(int Lid, int o, int n, int m, int M, double h, double k,
            double *p, double *g0, double *g1, double *u0, double *u1, double &F){
	double lhs, rhs;

	syst_lhs(Lid, o, n, m, M, h, k, p, g1, u1, lhs);	
	syst_rhs(Lid, o, n, m, M, h, k, p, g0, u0, rhs);	
	F   = lhs - rhs;
}

// parabolic operator F(u) at for all xm at fixed tn
void syst_F(int Lid, int o, int n, int M, double h, double k,
            double *p, double *g0, double *g1, double *u0, double *u1, double *F){
	int m;
	double Fm;
	for (m = 0; m < M; m++){
		syst_F(Lid, o, n, m, M, h, k, p, g0, g1, u0, u1, Fm);	
		F[m] = Fm;
	}
}

// left-hand side (unknowns) at (tn,xm)
void syst_lhs(int Lid, int o, int n, int m, int M, double h, double k,
              double *p, double *g1, double *u1, double &lhs){
	double cf = 0.5*k;	// trapezoidal weight
	double Lu1;					// spatial difference operator
	char   LR = 'L';
	
	syst_L(LR, Lid, o, n+1, m, M, h, k, p, u1, Lu1);
	lhs = u1[m] - cf*(Lu1 + g1[m]);
}

// right-hand side operator (knowns) at (tn,xm)
void syst_rhs(int Lid, int o, int n, int m, int M, double h, double k,
              double *p, double *g0, double *u0, double &rhs){
	double cf = 0.5*k;	// trapezoidal weight
	double Lu0;					// spatial difference operator
	char   LR = 'R';
	
	syst_L(LR, Lid, o, n  , m, M, h, k, p, u0, Lu0);
	rhs = u0[m] + cf*(Lu0 + g0[m]);
}

// partial derivative dF_i/du1_j
void syst_DF(int Lid, int o, int n, int mi, int mj, int M, double h, double k,
             double *p, double *g1, double *u1, double *lhs, double &DF){
	int m;
	double u1p[M];
	double lhsp  ;
	double du1 = 0.0001;

	for (m = 0; m < M; m++)
		u1p[m] = u1[m];
	u1p[mj] += du1;

	// calculate partial derivative by forward differences
	syst_lhs(Lid, o, n, mi, M, h, k, p, g1, u1p, lhsp);	
	DF = (lhsp - lhs[mi])/du1;
}

// Jacobian matrix DF = dF/du1 (all xmi, xmj)
//  NOTE: sprs = indicator function
void syst_DF(bool sprs, int Lid, int o, int n, int M, double h, double k,
             double *p, double *g1, double *u1, double *lhs, double *DF){
	int mi, mj;
	double DFij;
		
	if (sprs == false){
		/*------ CALCULATE PARTIAL DERIVATIVES ------*/
		/*------  WITHOUT EXPLOITING SPARSITY  ------*/
		for (mi = 0; mi < M; mi++){			// element of F
			for (mj = 0; mj < M; mj++){		// element of u
				syst_DF(Lid, o, n, mi, mj, M, h, k, p, g1, u1, lhs, DFij);
				DF[mi*M + mj] = DFij;
			}
		}
	}
	
	if (sprs == true){
		/*------ CALCULATE PARTIAL DERIVATIVES ------*/
		/*------      EXPLOITING SPARSITY      ------*/
		for (mi = 0; mi < M; mi++){
			for (mj = 0; mj < M; mj++){
				if (  ((o == 2) &&	// second-order spatial difference operator
					    (  (mi == 0   && mj < mi+4)		// forward differences
					    || (mi == M-1 && mj > mi-4)		// backward differences
					    || (mi > 0    && mi < M-1 		// central differences
					 	     && mj < mi+2 && mj > mi-2)
					 	  ))
					 || ((o == 3) &&	// third-order
					    (  (mi == 0   && mj < mi+5) 	// forward  differences
					    || (mi == M-1 && mj > mi-5)		// backward differences
					    || (mi > 0    && mi < M-1 
					 	 && mj < mi+2 && mj > mi-2)			// central differences
					 	  ))
					 || ((o == 4) &&	// fourth-order
				      (  (mi == 0   && mj < mi+6)		// forward differences
					    || (mi == 1   && mj < mi+6)
					    || (mi == M-1 && mj > mi-6)		// backward differences
					    || (mi == M-2 && mj > mi-6)
					    || (mi > 1    && mi < M-2 		// central differences
					 	     && mj < mi+3 && mj > mi-3)
					 	  ))
					)
						syst_DF(Lid, o, n, mi, mj, M, h, k, p, g1, u1, lhs, DFij);
					else
						DFij = 0.0;
				DF[mi*M + mj] = DFij;
			}
		}
	}
}

// spatial difference operator L(u) at (tn,xm)
//  NOTE: this is the only function that depends on Lid (problem type)
void syst_L(char LR, int Lid, int o, int n, int m, int M, double h, double k, 
						double *p, double *u, double &Lu){
	double cf           ; // coefficient
	double d1u          ;	// first -order spatial difference operators
	double d2u, l2u, b2u;	// second-     
	double d3u, l3u, b3u;	// third -     
	double d4u, l4u, b4u;	// fourth-     

	/* Heat equation in one dimension with constant diffusivity
	 * with Dirichlet boundary conditions. */
	if (Lid == 0){
		// difference operators
		diff_d2(m, M, h, u, d2u);

		// parameters
		cf = p[0];						// diffusion coefficient (constant)

		// heat equation operator
		Lu = d2u;

		if (LR == 'L' && (m == 0 || m == M-1))
			Lu = u[m];
		if (LR == 'R' && (m == 0 || m == M-1))
			Lu = 0.0 ;
	}

	/* Axisymmetric spreading of a thin film over a horizontal,
	 * planar substrate due to gravity. */
	if (Lid == 1){
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
	
		// STILL NEED BOUNDARY CONDITIONS HERE 
		// STILL NEED BOUNDARY CONDITIONS HERE 
		// STILL NEED BOUNDARY CONDITIONS HERE 
		// STILL NEED BOUNDARY CONDITIONS HERE 
		// STILL NEED BOUNDARY CONDITIONS HERE 

	}
}

#endif
