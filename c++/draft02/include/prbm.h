/* PROBLEM TYPE
 *  General parabolic equation:
 *          du
 *   F(u) = -- - L(u) - g = 0
 *          dt
 *  where L is an operator in one space dimension and g is a source term.
 *  In general, u = u(t,x;p) where t is time, x is the spatial coordinate,
 *  and p is a vector of parameters.
 *  Discretize time domain using the trapezoidal rule.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 *  u   [input]			unknown function
 *  L		[output]		spatial difference operator
 *  g   [output]		source function
 */

#ifndef PRBM_H
#define PRBM_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_bessel.h>
#include "./diff.h"

using namespace std;

/* PROTOTYPES */
void prbm_L (char, int, int, int, int, int, int, int, double, double, double *, double *, double &);
void prbm_L2(char, double, double, double, double, double &);
void prbm_BC(char, char, double, double, double &);

/* IMPLEMENTATIONS */

// compute L[u(t,x)] + g(t,x) at (tn,xm)
void prbm_L(char LR, int Lid, int dim, int ord, int var, int n, int m, int M, double h, double k, 
						double *p, double *u, double &Lug){
	int    i, j;
	double cf, cf2, cf3 ; // coefficients
	double d1u          ;	// first -order spatial difference operators
	double d2u, l2u, b2u;	// second-     
	double d3u, l3u, b3u;	// third -     
	double d4u, l4u, b4u;	// fourth-     
	
	double Lu, g        ;
	double x = m*h      ;
	double t = n*k      ;

	char   BC           ;	// BC type

	// difference operators
	diff_d1(m, M, h, u, d1u);
	diff_d2(m, M, h, u, d2u);
	diff_d3(m, M, h, u, d3u);
	diff_d4(m, M, h, u, d4u);

	diff_l2(m, M, h, u, l2u);

	diff_b2(m, M, h, u, b2u);
	diff_b3(m, M, h, u, b3u);
	diff_b4(m, M, h, u, b4u);
			
	// trapezoidal weight
	double w = 0.5*k;
	
	// parameters
	cf  = p[0];		// diffusivity (constant)
	cf2 = p[1];		// curvature of lower boundary
			
	// initialize source
	g = 0.0;

	/* SECOND-ORDER SPATIAL OPERATORS
	 *	Lid = 0 : diffusion (Laplacian operator)
	 *  Lid = 1 : gravitational spreading
	 */
	if (ord == 2){
		if (dim == 1){			// 1-D on a line
			if (Lid == 0){
				g  = 0.0;
				Lu = cf*d2u;
				prbm_L2(LR, w, g, u[m], Lu, Lug);
				if (m == 0)
					prbm_BC(LR, 'N', u[m], d1u, Lug);
				if (m == M-1)
					prbm_BC(LR, 'N', u[m], d1u, Lug);
			}
			if (Lid == 1){	
				cout << "Choose different problem" << endl;
				return;
			}
		}
		if (dim == 2){			// 1-D w/azimuthal symmetry about x = 0
			if (Lid == 0){
				g  = 0.0;
				Lu = cf*l2u;
				prbm_L2(LR, w, g, u[m], Lu, Lug);
				if (m == 0)
					prbm_BC(LR, 'N', u[m], d1u, Lug);
				if (m == M-1)
					prbm_BC(LR, 'N', u[m], d1u, Lug);
			}
			if (Lid == 1){
				cout << "Choose different problem" << endl;
				return;
			}
		}
	}

	/* FOURTH-ORDER SPATIAL OPERATORS
	 *  --- need at least 2 variables to properly distribute the BCs ---
	 *  Lid = 0 : capillary spreading
	 *  Lid = 1 : gravity-capillary spreading
	 *            [cf. Moriarty, Schwartz, and Tuck (1991)]
	 */
	if (ord == 4){
		if (dim == 1){			// 1-D on a line
			if (Lid == 0){
				cout << "Choose different problem" << endl;
				return;
			}
			if (Lid == 	1){
				g = 0.0;
				
				// assume M is even
				int m2;								// local index
				int M2 = M/2;					// upper bound for local index
				double H[M2], Q[M2];
				double d1H, d1Q;
				for (i = 0; i < M2; i++){
					H[i] = u[i   ];
					Q[i] = u[i+M2];
				}
				if (m < M2)			 {	// mass balance for H(x,t)
					if      (m == 0 )  {	// Dirichlet condition on h at m2 = 0
						diff_d1(0, M2, h, H, d1H);
						prbm_BC(LR, 'D', H[0], d1H, Lug);
					}
					else if (m == M2-1){	// Dirichlet condition on h at m2 = M2-1
						diff_d1(M2-1, M2, h, H, d1H);
						prbm_BC(LR, 'D', H[M2-1], d1H, Lug);
						if (LR == 'R')
							Lug = 0.1; // far-field thickness
					}
					else               {
						m2 = m;
						diff_d1(m2, M2, h, Q, d1Q);
						Lu = -d1Q;
						prbm_L2(LR, w, g, H[m2], Lu, Lug);
					}
				}
				else             {	// flux Q(x,t) (eq of motion + normal stress condition)
					if      (m == M2  ){				// Neumann condition on h at m2 = M2-1
						diff_d1(M2-1, M2, h, H, d1H);
						prbm_BC(LR, 'N', H[M2-1], d1H, Lug);
					}
					else if (m == M -1){				// Dirichlet  condition on q at m2 = 0
						diff_d1(0, M2, h, Q, d1Q);
						prbm_BC(LR, 'D', Q[0], d1Q, Lug); 
					}
					else               {
						if      (LR == 'L'){
							double H2, H3, d3H;
							m2 = m - M2    ; 	// subtract 1 b/c we placed both BCs at top
							H2 = H[m2]*H[m2];
							H3 = H[m2]*H3  ;
							diff_d3(m2, M2, h, H, d3H);
							Lu = (1.0/3.0)*H3*(1.0 + cf*d3H);
							Lug = Q[m2] - Lu;
						}
						else if (LR == 'R')
							Lug = 0.0;
					}
				}
			}
		}
		if (dim == 2){			// 1-D w/azimuthal symmetry about x = 0
			if (Lid == 0){
				cout << "Choose different problem" << endl;
				return;
			}
			if (Lid == 1){
				cout << "Choose different problem" << endl;
				return;
			}
		}
	}
}

void prbm_L2(char LR, double w, double g, double u, double Lu, double &Lug){
	if      (LR == 'L')
		Lug = u - w*(Lu + g);
	else if (LR == 'R')
		Lug = u + w*(Lu + g);
}

void prbm_BC(char LR, char BC, double u, double d1u, double &Lug){
	if      (LR == 'L'){
		if      (BC == 'D'){
			Lug = u;
		}
		else if (BC == 'N'){
			Lug = d1u;
		}
	}
	else if (LR == 'R')
		Lug = 0.0;
}

#endif
