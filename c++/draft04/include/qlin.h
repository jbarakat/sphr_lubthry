/* QUASI-LINEARIZATION
 *  Solve the quasi-linearized equations by a predictor-corrector method.
 *
 * REFERENCES
 *  von Rosenberg, Methods for the Numerical Solution of PDEs (1969)
 *  Lopez et al, J Coll Int Sci (1976) - gravitational spreading
 *  Moriarty, Schwartz, and Tuck, Phys Fluids A (1991) - gravity-draining w/surface tension
 *  
 * PARAMETERS
 *  n    [input]   time index that ranges from 0 to N-1
 *  m    [input]   space index that ranges from 0 to M-1
 *  N    [input]   number of time steps
 *  M    [input]   number of grid points
 *  k    [input]   time spacing
 *  h    [input]   grid spacing
 *  var  [input]   number of variables
 *  u    [input]   unknown function (size = M*var)
 *  p    [input]   parameters (arbitrary size)
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
#include "./diff.h"

using namespace std;

/* PROTOTYPES */

/* IMPLEMENTATIONS */

// Project u0 at t = t_n to t = t_{n+1/2}
void qlin_proj(int n, int M, double dx, double dt, double *p, double *u0, double *u1){
	int var = 2;

	// system size
	int S = var*M;



}

// Correct solution at t = t_{n+1}
void qlin_corr(){
}

// calculate F and DF
void syst_F_DF(int sprs, int dim, int var, int n, int M, double h, double k,
	             double *p, double *u0, double *u1, double *F, double *DF){
	syst_F (      dim, var, n, M, h, k, p, u0, u1, F);
	syst_DF(sprs, dim, var, n, M, h, k, p, u0, u1, F, DF);
}

// calculate F
void syst_F(int dim, int var, int n, int m, int M, double h, double k,
            double *p, double *u0, double *u1, double *F){
	// ASSUME HAVE 2 VARIABLES
	if (var == 2){
		// arrange dependent variables as such:
		//  u = [H0, Q0, H1, Q1, ..., Hm, Qm, ..., HM-1, QM-1]
		int i;
		double H1[M], Q1[M];
		double H0[M], Q0[M];
		double lhsH, rhsH, lhsQ, rhsQ;

		// initialize H, Q
		for (i = 0; i < M; i++){
			H0[i] = u0[2*i  ];
			H1[i] = u1[2*i  ];
			Q0[i] = u0[2*i+1];
			Q1[i] = u1[2*i+1];
		}
		
		// boundary conditions
		if      (m == 0  ){	
			lhsH = H1[m]; rhsH = 0.0;
			lhsQ = Q1[m]; rhsQ = 0.0;
		}
		else if (m == M-1){
			double d1H1;
			double cf2 = p[1];	// far-field thickness

			diff_d1(m, M, h, H1, d1H1);
			
			lhsH = H1[m]; rhsH = cf2;
			lhsQ = d1H1 ; rhsQ = 0.0;
		}
		
		// interior nodes
		else              {
			double d1Q0, d1Q1, d3H1, H12, H13;
			double w = 0.5*k; 		// trapezoidal weight
			double cf1 = p[0];		// reciprocal Bond number

			diff_d1(m, M, h, Q0, d1Q0);
			diff_d1(m, M, h, Q1, d1Q1);
			diff_d3(m, M, h, H1, d3H1);
			H12 = H1[m]*H1[m];
			H13 = H12  *H1[m];
			
			// mass balance
			lhsH = H1[m] + w*d1Q1;
			rhsH = H0[m] - w*d1Q0;
			
			// momentum balance + normal stress balance
			lhsQ = Q1[m] - (1.0/3.0)*H13*(1.0 + cf1*d3H1);
			rhsQ = 0.0;
		}

		// assemble (nonlinear) difference equation
		F[0] = lhsH - rhsH;
		F[1] = lhsQ - rhsQ;
	}
	else {
		cout << "Should have 2 variables H and Q." << endl;
	}
}

// calculate F
void syst_F(int dim, int var, int n, int M, double h, double k,
            double *p, double *u0, double *u1, double *F){

	int m, v;
	double Fm[var];
	
	for (m = 0; m < M; m++){
		syst_F(dim, var, n, m, M, h, k, p, u0, u1, Fm);
		for (v = 0; v < var; v++){
			F[var*m + v] = Fm[v];
		}
	}
}

// Jacobian matrix DF = dF/du1 (all xmi, xmj)
//  NOTE: sprs = sparsity pattern
void syst_DF(int sprs, int dim, int var, int n, int M, double h, double k,
             double *p, double *u0, double *u1, double *F, double *DF){
	int m, v, i, j;
	int Mv = var*M;
	double u1p[Mv];
	double Fp[var], DFp[var];
	double du1 = 0.0001;

	/*------ CALCULATE PARTIAL DERIVATIVES ------*/
	/*------  WITHOUT EXPLOITING SPARSITY  ------*/
	if (sprs == 0){
		for (m = 0; m < M; m++){			// element of F
			for (i = 0; i < Mv; i++){		// element of u
				// perturb u1
				for (j = 0; j < Mv; j++)
					u1p[j] = u1[j];
				u1p[i] += du1;

				syst_F(dim, var, n, m, M, h, k, p, u0, u1p, Fp);
				for (v = 0; v < var; v++){
					DFp[v]                  = (Fp[v] - F[var*m + v])/du1;
					DF [(var*m + v)*Mv + i] = DFp[v];
				}
			}
		}
	}
	
	/*------ CALCULATE PARTIAL DERIVATIVES ------*/
	/*------      EXPLOITING SPARSITY      ------*/
}

#endif
