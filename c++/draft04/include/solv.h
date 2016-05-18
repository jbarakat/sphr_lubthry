/* NUMERICAL SOLVER
 *  Assemble the difference system using the trapezoidal method (Crank-Nicholson scheme).
 *
 * REFERENCES
 *  Lopez et al, J Coll Int Sci (1976) - gravitational spreading
 *  Moriarty, Tuck, and Schwartz, Phys Fluids A (1991) - gravitational thinning w/surface tension
 *  
 * PARAMETERS
 *  M		[input]			number of grid points
 *  p		[input]			parameters
 *  u		[input]			solution vector
 */

#ifndef SOLV_H
#define SOLV_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include "./diff.h"
#include "./lalg.h"

using namespace std;

/* PROTOTYPES */

/* IMPLEMENTATIONS */

// prediction
void solv_pred(int M, double dx, double dt, double *p, double *u0, double *u1){
	int    i, j;
	int    R = M-2; // system size
	double a0[R], b0[R], c0[R], d0[R], e0[R], f0[R];
	double a [R], b [R], c [R], d [R], e [R], f [R];
	double k0[R], k[R];
	double v0[R], v[R];
	double dv;
	double tol = 1e-6;

	// parameters
	double B   = p[0];	// Bond number
	double U   = p[1];	// far-field thickness
	double cf1 = 0.25*dt/dx;
	double cf2 = cf1/(B*dx*dx*dx);

	// initialize
	for (i = 0; i < R; i++){
		v0[i] = u0[i];
		k0[i] = pow(0.5*(u0[i] + u0[i+1]),3)/3.0;
	}
	k0[R-1] = pow(0.5*(u0[i] + U),3)/3.0;

	a0[0] = 0.0;
	b0[0] = 0.0;
	c0[0] = 
	d0[0] = 
	e0[0] = 
	f0[0] = 

	a0[1] = 0.0;
	b0[1] = 
	c0[1] = 
	d0[1] = 
	e0[1] = 
	f0[1] = 

	for (i = 2; i < R-2; i++){
		a0[i] =       cf2*     k0[i-1]             ;
		b0[i] =     - cf2*(3.0*k0[i-1] +     k0[i]);
		c0[i] = 1.0 - cf2*(3.0*k0[i-1] + 3.0*k0[i]);
		d0[i] =     - cf2*(    k0[i-1] + 3.0*k0[i]);
		e0[i] =       cf2*                   k0[i] ;
		f0[i] = 
	}

	a0[R-2] = 
	b0[R-2] = 
	c0[R-2] = 
	d0[R-2] = 
	e0[R-2] = 0.0;
	f0[R-2] = 
	
	a0[R-1] = 
	b0[R-1] = 
	c0[R-1] = 
	d0[R-1] = 0.0;
	e0[R-1] = 0.0;
	f0[R-1] = 
	
	// iterate
	for (i = 0; i < R-1; i++)
		k[i] = pow(0.25*(u0[i  ] + v0[i  ] + u0[i+1] + v0[i+1]),3)/3.0;
	k[R-1] = pow(0.25*(u0[R-1] + v0[R-1] + 2.0*U),3)/3.0;

	// compute coefficients for pentadiagonal system
	a[0]   = 0.0               ;
	b[0]   = 0.0               ;
	c[0]   = 1.0 + 4.0*cf2*k[0];
	d[0]   =     - 3.0*cf2*k[0];
	e[0]   =           cf2*k[0]; 
	f[0]   = 

	a[1]   = 0.0                            ;
	b[1]   =     - cf2*(3.0*k[0] +     k[1]);
	c[1]   = 1.0 + cf2*(3.0*k[0] + 3.0*k[1]);
	d[1]   =     - cf2*(    k[0] + 3.0*k[1]);
	e[1]   =       cf2*                k[1] ;
	f[1]   = 

	for (i = 2; i < R-2; i++){
		a[i] =       cf2*     k[i-1]            ;
		b[i] =     - cf2*(3.0*k[i-1] +     k[i]);
		c[i] = 1.0 + cf2*(3.0*k[i-1] + 3.0*k[i]);
		d[i] =     - cf2*(    k[i-1] + 3.0*k[i]);
		e[i] =       cf2*                  k[i] ;
		f[i] = 
	}

	a[R-2] =       cf2*     k[R-3]              ;
	b[R-2] =     - cf2*(3.0*k[R-3] +     k[R-2]);
	c[R-2] = 1.0 + cf2*(3.0*k[R-3] + 3.0*k[R-2]);
	d[R-2] =     - cf2*(    k[R-3] + 3.0*k[R-2]);
	e[R-2] = 0.0                                ;
	f[R-2] = 
	
	a[R-1] =       cf2*     k[R-2]              ;
	b[R-1] =     - cf2*(3.0*k[R-2] +     k[R-1]);
	c[R-1] = 1.0 + cf2*(3.0*k[R-2] + 3.0*k[R-1]);
	d[R-1] = 0.0                                ;
	e[R-1] = 0.0                                ;
	f[R-1] = 

	// solve the pentadiagonal system for v

	// compute norm of difference ||v - v0||
	
}

// correction
void solv_corr(){
}


#endif
