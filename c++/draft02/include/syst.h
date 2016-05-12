/* DIFFERENCE SYSTEM
 *  Assemble the difference system using the trapezoidal method (Crank
 *  -Nicholson scheme).
 *
 * REFERENCES
 *  Lopez et al, J Coll Int Sci (1976) - gravitational spreading
 *  
 * PARAMETERS
 *  n    [input]   time index that ranges from 0 to N-1
 *  m    [input]   space index that ranges from 0 to M-1
 *  N    [input]   number of time steps
 *  M    [input]   number of grid points (system size)
 *  k    [input]   time spacing
 *  h    [input]   grid spacing
 *  u    [input]   unknown function (size = M)
 *  g    [input]   (known) source function (size = M)
 *  p    [input]   parameters (arbitrary size)
 *  ord  [input]   order of the spatial difference operator L
 *  Lid  [input]	 indicator for problem type
 *  sprs [input]	 indicator for whether or not to exploit sparsity
 *  F		 [output]	 function to be set to zero
 *  DF	 [output]	 Jacobian matrix
 */

#ifndef SYST_H
#define SYST_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include "./prbm.h"

/* PROTOTYPES */
void syst_F_DF(bool, int, int, int, int, int,           int, double, double, double *, double *, double *, double *, double *);
void syst_F   (      int, int, int, int, int, int,      int, double, double, double *, double *, double *, double &          );
void syst_F   (      int, int, int, int, int,           int, double, double, double *, double *, double *, double *          );
void syst_DF  (      int, int, int, int, int, int, int, int, double, double, double *,           double *, double *, double &);
void syst_DF  (bool, int, int, int, int, int,           int, double, double, double *,           double *, double *, double *);
void syst_lhs (      int, int, int, int, int, int,      int, double, double, double *,           double *, double &          );
void syst_rhs (      int, int, int, int, int, int,      int, double, double, double *, double *,           double &          );

/* IMPLEMENTATIONS */

// assemble F and DF for all xm at fixed tn
void syst_F_DF(bool sprs, int Lid, int dim, int ord, int var, int n, int M, double h, double k,
	             double *p, double *u0, double *u1, double *F, double *DF){
	int m;
	double Fm, lhsm, rhsm;
	double lhs[M], rhs[M], DFm[M*M];

	// calculate F
	for (m = 0; m < M; m++){
		syst_lhs(Lid, dim, ord, var, n, m, M, h, k, p, u1, lhsm);	
		syst_rhs(Lid, dim, ord, var, n, m, M, h, k, p, u0, rhsm);	
		Fm   = lhsm - rhsm;

		lhs[m] = lhsm;
		rhs[m] = rhsm;
		F  [m] = Fm  ;
	}
	
	// calculate DF
	syst_DF(sprs, Lid, dim, ord, var, n, M, h, k, p, u1, lhs, DFm);
	for (m = 0; m < M*M; m++){
		DF[m] = DFm[m];
	}
}

// parabolic operator F(u) at (tn,xm)
// time difference using the trapezoidal method (Crank-Nicholson)
void syst_F(int Lid, int dim, int ord, int var, int n, int m, int M, double h, double k,
            double *p, double *u0, double *u1, double &F){
	double lhs, rhs;

	syst_lhs(Lid, dim, ord, var, n, m, M, h, k, p, u1, lhs);	
	syst_rhs(Lid, dim, ord, var, n, m, M, h, k, p, u0, rhs);	
	F   = lhs - rhs;
}

// parabolic operator F(u) at for all xm at fixed tn
void syst_F(int Lid, int dim, int ord, int var, int n, int M, double h, double k,
            double *p, double *u0, double *u1, double *F){
	int m;
	double Fm;
	for (m = 0; m < M; m++){
		syst_F(Lid, dim, ord, var, n, m, M, h, k, p, u0, u1, Fm);	
		F[m] = Fm;
	}
}

// partial derivative dF_i/du1_j
void syst_DF(int Lid, int dim, int ord, int var, int n, int mi, int mj, int M, double h, double k,
             double *p, double *u1, double *lhs, double &DF){
	int m;
	double u1p[M];
	double lhsp  ;
	double du1 = 0.0001;

	for (m = 0; m < M; m++)
		u1p[m] = u1[m];
	u1p[mj] += du1;

	// calculate partial derivative by forward differences
	syst_lhs(Lid, dim, ord, var, n, mi, M, h, k, p, u1p, lhsp);	
	DF = (lhsp - lhs[mi])/du1;
}

// Jacobian matrix DF = dF/du1 (all xmi, xmj)
//  NOTE: sprs = indicator function
void syst_DF(bool sprs, int Lid, int dim, int ord, int var, int n, int M, double h, double k,
             double *p, double *u1, double *lhs, double *DF){
	int mi, mj;
	double DFij;
		
	if (sprs == false){
		/*------ CALCULATE PARTIAL DERIVATIVES ------*/
		/*------  WITHOUT EXPLOITING SPARSITY  ------*/
		for (mi = 0; mi < M; mi++){			// element of F
			for (mj = 0; mj < M; mj++){		// element of u
				syst_DF(Lid, dim, ord, var, n, mi, mj, M, h, k, p, u1, lhs, DFij);
				DF[mi*M + mj] = DFij;
			}
		}
	}
	
	if (sprs == true){
		/*------ CALCULATE PARTIAL DERIVATIVES ------*/
		/*------      EXPLOITING SPARSITY      ------*/
		for (mi = 0; mi < M; mi++){
			for (mj = 0; mj < M; mj++){
				if (  ((ord == 2) &&	// second-order spatial difference operator
					    (  (mi == 0   && mj < mi+4)		// forward differences
					    || (mi == M-1 && mj > mi-4)		// backward differences
					    || (mi > 0    && mi < M-1 		// central differences
					 	     && mj < mi+2 && mj > mi-2)
					 	  ))
					 || ((ord == 3) &&	// third-order
					    (  (mi == 0   && mj < mi+5) 	// forward  differences
					    || (mi == M-1 && mj > mi-5)		// backward differences
					    || (mi > 0    && mi < M-1 
					 	 && mj < mi+2 && mj > mi-2)			// central differences
					 	  ))
					 || ((ord == 4) &&	// fourth-order
				      (  (mi == 0   && mj < mi+6)		// forward differences
					    || (mi == 1   && mj < mi+6)
					    || (mi == M-1 && mj > mi-6)		// backward differences
					    || (mi == M-2 && mj > mi-6)
					    || (mi > 1    && mi < M-2 		// central differences
					 	     && mj < mi+3 && mj > mi-3)
					 	  ))
					)
						syst_DF(Lid, dim, ord, var, n, mi, mj, M, h, k, p, u1, lhs, DFij);
					else
						DFij = 0.0;
				DF[mi*M + mj] = DFij;
			}
		}
	}
}

// left-hand side (unknowns) at (tn,xm)
void syst_lhs(int Lid, int dim, int ord, int var, int n, int m, int M, double h, double k,
              double *p, double *u1, double &lhs){
	double Lu1;					// = u1 - L(u1) + g1
	char   LR = 'L';
	
	prbm_L(LR, Lid, dim, ord, var, n+1, m, M, h, k, p, u1, Lu1);
	lhs = Lu1;
}

// right-hand side operator (knowns) at (tn,xm)
void syst_rhs(int Lid, int dim, int ord, int var, int n, int m, int M, double h, double k,
              double *p, double *u0, double &rhs){
	double Lu0;					// = u0 + L(u0) + g0
	char   LR = 'R';
	
	prbm_L(LR, Lid, dim, ord, var, n  , m, M, h, k, p, u0, Lu0);
	rhs = Lu0;
}

#endif
