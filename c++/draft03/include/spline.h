/* INTERPOLATION
 *  Generate interpolants using cubic splines.
 *   x(s) - x0 = a*(s - s0)^3 + b*(s - s0)^2 + c*(s - s0)
 *   y(s) - y0 = a*(s - s0)^3 + b*(s - s0)^2 + c*(s - s0)
 *
 * REFERENCES
 *  Moin, Fundamentals of Engineering Numerical Analysis, Cambridge University Press (2010) (Ch. 1)
 *  Pozrikidis, A Practical Guide to Boundary Element Methods, Chapman & Hall/CRC (2002) (Ch. 3)
 *  
 * PARAMETERS
 *  x,y    [input]		set of N+1 grid points
 *  xi,yi  [output] 	interpolated grid point
 *  N      [input]		number of grid points
 *  a      [output]		spline coefficient of 3rd derivative
 *  b      [output]		spline coefficient of 2nd derivative
 *  c      [output]		spline coefficient of 1st derivative
 */

#ifndef SPLINE_H
#define SPLINE_H

/* HEADER FILES */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <lapacke.h>
#include <gsl/gsl_sf_trig.h>

/* PROTOTYPES */
void spline(const int, double*, double*,
            double, double, double, double,
            double*, double*, double*,
            double*, double*, double*);
void spline(const int, double*, double*, double, double,
            double*, double*, double*);
void spline(const int, double*, double*, double, double,
            double, double&);



/* IMPLEMENTATIONS */

/* Generate 3N coefficients a, b, c for cubic spline interpolation of
 * N+1 points (x,y(x)) by solving a tridiagonal system of equations with
 * clamped boundary conditions.
 */
void spline(const int N, double *x, double *y, double slope1, double slope2,
            double *a, double *b, double *c){
	// declare variables
	int i, j, info;
	double h [N  ], rhs[N-1], sln[N-1];
	double dl[N-2], d  [N-1], du [N-2];
	
	// compute intervals h
	for (i = 0; i < N; i++){
		h[i] = x[i+1] - x[i];
	}

	// generate tridiagonal matrix for N-1 equations
	d [0] = 2.*(h[0] + h[1]) - 0.5*h[0];
	dl[0] = h[1];
	du[0] = h[1];

	for (i = 1; i < N-2; i++){
		d [i] = 2.*(h[i] + h[i+1]);
		dl[i] = h[i+1];
		du[i] = h[i+1];
	}

	d[N-2] = 2.*(h[N-2] + h[N-1]) - 0.5*h[N-1];

	for (i = 0; i < N-1; i++){
		rhs[i] = 3.*((y[i+2] - y[i+1])/h[i+1] - (y[i+1] - y[i])/h[i]);
	}

	rhs[0]   -= 1.5*((y[1] - y[0])  /h[0]   - slope1);
	rhs[N-2] += 1.5*((y[N] - y[N-1])/h[N-1] - slope2);

	// solve the tridiagonal system
	info = LAPACKE_dgtsv(LAPACK_COL_MAJOR, N-1, 1, dl, d, du, rhs, N-1);
	
	for (i = 0; i < N-1; i++){
		sln[i] = rhs[i];
	}

	// compute the coefficients a, b, c
	for (i = 0; i < N-1; i++){
		b[i+1] = sln[i];
	}

	b[0] = -0.5*b[1]   + 1.5*((y[1] - y[0]  )/h[0]   - slope1)/h[0];
//	b[N] = -0.5*b[N-1] - 1.5*((y[N] - y[N-1])/h[N-1] - slope2)/h[N-1];
	
	for (i = 0; i < N; i++){
		a[i] = (b[i+1] - b[i])/(3.*h[i]);
		c[i] = (y[i+1] - y[i])/h[i] - h[i]*(b[i+1] + 2.*b[i])/3.;
	}
}

/* Generate 6N coefficients a, b, c for cubic spline interpolation of
 * N+1 points (x(s),y(s)), parametrized by the polygonal arc length s.
 */
void spline(const int N, double *x, double *y,
            double slopex1, double slopex2,
            double slopey1, double slopey2,
						double *ax, double *bx, double *cx,
						double *ay, double *by, double *cy){
	// declare variables
  int i, j, n;
	double s[N+1];
	double dx, dy, ds;

	// compute the polygonal arc length
	s[0] = 0;
	for (i = 1; i < N+1; i++){
		dx = x[i] - x[i-1];
		dy = y[i] - y[i-1];
		ds = sqrt(dx*dx + dy*dy);
		s[i] = s[i-1] + ds;
	}

	// compute cubic spline coefficients
	spline(N, s, x, slopex1, slopex2, ax, bx, cx);
	spline(N, s, y, slopey1, slopey2, ay, by, cy);
}

/* Interpolate to point yi = y(xi) based on a set of N+1 nodes (x,y(x)) using
 * cubic splines.
 */
void spline(const int N, double *x, double *y, double slope1, double slope2,
            double xi, double &yi){
  // declare variables
  int i, j, n;
  double a[N], b[N], c[N];
  double dx;

  // compute cubic spline coefficients
  spline(N, x, y, slope1, slope2, a, b, c);

  // find where xi lies in the domain 
  n = -1;
  for (i = 0; i < N; i++){
    if (xi > x[i] && xi < x[i+1]){
      n = i;
      break;
    }
    
    if (xi < x[i] && xi > x[i+1]){
      n = i;
      break;
    }
  }
  
  if (n < 0){
    printf("Error: xi not in x.");
    return;
  }

  // compute yi
  dx = xi - x[n];
  yi = x[n] + ((a[n]*dx + b[n])*dx + c[n])*dx;
}


#endif
