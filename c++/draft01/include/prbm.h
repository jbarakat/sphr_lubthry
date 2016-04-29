/* PROBLEM TYPE
 *  General parabolic equation:
 *          du
 *   F(u) = -- - L(u) - g = 0
 *          dt
 *  where L is an operator in one space dimension and g is a source term.
 *  In general, u = u(t,x;p) where t is time, x is the spatial coordinate,
 *  and p is a vector of parameters.
 *
 * REFERENCES
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
void prbm_L(char, int, int, int, int, double, double, double *, double *, double &);
void prbm_g(      int, int, int, int, double, double, double *,           double &);

/* IMPLEMENTATIONS */

// spatial difference operator L(u) at (tn,xm)
void prbm_L(char LR, int Lid, int n, int m, int M, double h, double k, 
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

// source function g(tn,xm) 
void prbm_g(int gid, int n, int m, int M, double h, double k, double *p, double &g){
	if (gid == 0){ // constant source of unit strength
		g = 1.0;
	}
	if (gid == 1){ // decaying sinusoidal wave
		g = ((M_PI*M_PI) - 1.0)*gsl_sf_exp(-n*k)*gsl_sf_sin(M_PI*m*h);
	}
}

#endif
