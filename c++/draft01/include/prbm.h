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

/* PROTOTYPES */
void prbm_L(char, int, int, int, int, int, int, double, double, double *, double *, double &);

/* IMPLEMENTATIONS */

// compute L[u(t,x)] + g(t,x) at (tn,xm)
void prbm_L(char LR, int Lid, int D, int o, int n, int m, int M, double h, double k, 
						double *p, double *u, double &Lug){
	double cf, cf2, cf3 ; // coefficients
	double d1u          ;	// first -order spatial difference operators
	double d2u, l2u, b2u;	// second-     
	double d3u, l3u, b3u;	// third -     
	double d4u, l4u, b4u;	// fourth-     
	
	double Lu, g        ;
	double x = m*h      ;
	double t = n*k      ;

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
			
	// source
	g = ((M_PI*M_PI) - 1.0)*gsl_sf_exp(-t)*gsl_sf_sin(M_PI*x);
	g = 0.0;

	/* SECOND-ORDER SPATIAL OPERATORS
	 *	Lid = 0 : diffusion (Laplacian operator)
	 *  Lid = 1 : gravitational spreading
	 */
	if (o == 2){
		if (D == 1){			// 1-D on a line
			if (Lid == 0)
				Lu = cf*d2u;
//			if (Lid == 1)
		}
		if (D == 2){			// 1-D w/azimuthal symmetry about x = 0
			if (Lid == 0)
				Lu = cf*l2u;
			if (Lid == 1)
				Lu = cf*u[m]*u[m]*(u[m]*l2u/3.0 + d1u*d1u);
		}
	}

	/* FOURTH-ORDER SPATIAL OPERATORS
	 *  Lid = 0 : capillary spreading
	 *  Lid = 1 : gravity-capillary spreading
	 *  Lid = 2 : penetration of rigid sphere through initially flat interface
	 *  Lid = 3 : drainage of thin film over rigid sphere
	 */
	if (o == 4){
//		if (D == 1){
//		}
		if (D == 2){			// 1-D w/azimuthal symmetry about x = 0
//			if (Lid == 0)
//			if (Lid == 1)
			if (Lid == 2){
				// REVISE THIS LATER
				double us  = 1.0 - t + 0.5*cf2*x*x;	// approximation (parabola)
				double f   = u[m] + us;
				double d1f = d1u + cf2*x;						// approximation (parabola)
				double f2  = f*f;
				double f3  = f2*f;
				double V   = 1.0 + f3*b4u/3.0 + f2*b3u*d1f;		// axial velocity at free surface
				double U   = -f2*b3u/2.0;											// radial velocity at free surface

				Lu = V - U*d1u;																// kinematic condition
			}
		}
	}
			
	// interior nodes
	if (LR == 'L')
		Lug = u[m] - w*(Lu + g);
	if (LR == 'R')
		Lug = u[m] + w*(Lu + g);
	
	// boundary nodes
	// NEED TO PASS INDICATOR FUNCTION FOR TYPE OF BCs AT EACH BOUNDARY
	// NEED TO PASS INDICATOR FUNCTION FOR TYPE OF BCs AT EACH BOUNDARY
	// NEED TO PASS INDICATOR FUNCTION FOR TYPE OF BCs AT EACH BOUNDARY
	if (o > 1){
		if (LR == 'L' && (m == 0 || m == M-1)){
			if (m == 0){
				Lug = u[m];		// Dirichlet BC
				Lug = d1u;		// Neumann BC
			}
			if (m == M-1){
				Lug = u[m];		// Dirichlet BC
				Lug = d1u;		// Neumann BC
			}
		}
		if (LR == 'R' && (m == 0 || m == M-1))
			Lug = 0.0 ;
	}

	// REVISE THIS LATER 
	// REVISE THIS LATER 
	// REVISE THIS LATER 
	if (o == 4){ // REVISE LATER
		if (LR == 'L' && (m == 1 || m == M-2)){
			if (m == 1)
				Lug = d3u;	// symmetry requirement
			if (m == M-2){
				Lug = d1u;		// Neumann BC (NOT SURE IF THIS IS CORRECT)???
			}
		}
		if (LR == 'R' && (m == 1 || m == M-2))
			Lug = 0.0 ;
	}
}

#endif
