/* NUMERICAL SOLVER
 *  Assemble the difference system using the trapezoidal method (Crank-Nicholson scheme).
 *
 * REFERENCES
 *  Lopez et al, J Coll Int Sci (1976) - gravitational spreading
 *  Moriarty, Tuck, and Schwartz, Phys Fluids A (1991) - gravitational thinning w/surface tension
 *  von Rosenberg, Methods for the Numerical Solution of PDEs (1969)
 *  
 * PARAMETERS
 *  J		[input]			number of grid points
 *  p		[input]			parameters
 *  u		[input]			solution vector
 */

/* The following problem is adapted from Moriarty, Tuck, and Schwartz, Physics
 * of Fluids A (1991). The differential system is
 *   u_t = - q_x,
 *     q = (u^3/3)*(1 + (1/B)*u_xxx)
 * where u is the film thickness, q is the flux, B is the Bond number.
 * 
 * NOTE: The independent variable u is shifted 1/2 step from the nodal points,
 *       so that u_{j+1/2,n} = u(x_{j}, t_{n}). There are J+1 nodal points x_{j}
 *       where j = 0, 1, ..., J-1, J.
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
void solv_tstep (int, double, double, double *, double *, double *);
void solv_coeffs(int, double, double, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *);

/* IMPLEMENTATIONS */

// prediction
void solv_tstep(int J, double dt, double dx, double *p, double *u0, double *u1){
	int i, j;
	int R = J - 2; // system size
	double v0[R], v1[R];
	double a[R], b[R], c[R], d[R], e[R], f[R];

	// initialize
	for (j = 1; j < J-1; j++){
		i = j-1;
		v0[i] = u0[j];
	}
	
	// prediction: project v0 onto intermediate time point to get v1
	solv_coeffs(J, 0.5*dt, dx, p, v0, v0, a, b, c, d, e, f);
	lalg_pent  (R, a, b, c, d, e, f, v1);

	// correction: correct v1 by taking the full time step
	solv_coeffs(J,     dt, dx, p, v0, v1, a, b, c, d, e, f);
	lalg_pent  (R, a, b, c, d, e, f, v1)
	
	// copy solution
	for (j = 1; j < J-1; j++){
		i = j-1;
		u1[j] = v0[i];
	}
	u[0  ] = U0;
	u[J-1] = U1;
	u[J  ] = U1;	// extra value b/c function was shifted 1/2 step
}

// compute pentadiagonal coefficients
void solv_coeffs(int J, double dt, double dx, double *p, double *u0, double *v,
                 double *a, double *b, double *c, double *d, double *e, double *f){
	int i, j;
	int R = J - 2; // system size

	// parameters
	double B  = p[0];		// Bond number
	double U0 = p[1];		// left boundary value
	double U1 = p[2];		// right boundary value
	
	// setup
	double cf1 = dt/dx;
	double cf2 = 1.0/(B*dx*dx*dx);
	double g[R], g1, g2;
	double A, B, C, D, E, F;
	for (i = 0; i < R-1; i++)
		g[i] = (cf1/3.0)*pow((v[i] + v[i+1])/2.0,3);
	g[R-1] = (cf1/3.0)*pow((v[R-1] + U1)/2.0,3);

	// compute coefficients for the particular problem
	/*---------------*/
	j = 1  ; i = j-1;
	g1   = 0.0;
	g2   = cf2*g[i  ];
	
	C    =            4.0*g2 ;
	D    = -          3.0*g2 ;
	E    =                g2 ;
	F    = g[i]              ;

	a[i] = 0.0;
	b[i] = 0.0;
	c[i] = 1.0 + C;
	d[i] = D;
	e[i] = E;
	f[i] = (1.0 - C)*u0[j] - D*u0[j+1] - E*u0[j+2] - F;

	/*---------------*/
	j = 2  ; i = j-1;
	g1 = cf2*g[i-1];
	g2 = cf2*g[i  ];

	B    = -(4.0*g1 +     g2);
	C    =  (3.0*g1 + 3.0*g2);
	D    = -(    g1 + 3.0*g2);
	E    =                g2 ;
	F    = g[i] - g[i-1]     ;

	a[i] = 0.0;
	b[i] = B;
	c[i] = 1.0 + C;
	d[i] = D;
	e[i] = E;
	f[i] = -B*u0[j-1] + (1.0 - C)*u0[j] - D*u0[j+1] - E*u0[j+2] - F;

	/*---------------*/
	for (j = 3; j < J-3; j++){
		i = j-1;
		g1 = cf2*g[i-1];
		g2 = cf2*g[i  ];

		A    =       g1          ;
		B    = -(3.0*g1 +     g2);
		C    =  (3.0*g1 + 3.0*g2);
		D    = -(    g1 + 3.0*g2);
		E    =                g2 ;
		F    = g[i] - g[i-1]     ;

		a[i] = A;
		b[i] = B;
		c[i] = 1.0 + C;
		d[i] = D;
		e[i] = E;
		f[i] = -A*u0[j-2] - B*u0[j-1] + (1.0 - C)*u0[j] - D*u0[j+1] - E*u0[j+2] - F
	}
	
	/*---------------*/
	j = J-2; i = j-1;
	g1 = cf2*g[i-1];
	g2 = cf2*g[i  ];

	A    =       g1          ;
	B    = -(3.0*g1 +     g2);
	C    =  (3.0*g1 + 3.0*g2);
	D    = -(    g1 + 3.0*g2);
	F    = g[i] - g[i-1] + 2.0*U1;

	a[i] = A;
	b[i] = B;
	c[i] = 1.0 + C;
	d[i] = D;
	e[i] = 0.0;
	f[i] = -A*u0[j-2] - B*u0[j-1] + (1.0 - C)*u0[j] - D*u0[j+1] - F;

	/*---------------*/
	j = J-1; i = j-1;
	g1 = cf2*g[i-1];
	g2 = cf2*g[i  ];

	A    =       g1          ;
	B    = -(3.0*g1 +     g2);
	C    =  (3.0*g1 + 3.0*g2);
	F    = g[i] - g[i-1] + 2.0*U1 - 2.0*(g1 + 2.0*g2);

	a[i] = A;
	b[i] = B;
	c[i] = 1.0 + C;
	d[i] = 0.0;
	e[i] = 0.0;
	f[i] = -A*u0[j-2] - B*u0[j-1] + (1.0 - C)*u0[j] - F;
}


#endif
