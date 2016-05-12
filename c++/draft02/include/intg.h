/* TIME INTEGRATION
 *
 * REFERENCES
 *  
 * PARAMETERS
 */

#ifndef INTG_H
#define INTG_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include "./newt.h"

/* PROTOTYPES */
void intg_timeStep(bool, int, int, int, int, int, int, double, double, double *, double *, double *);
void intg_timeEvol(bool, int, int, int, int, int, int, double, double, double *, double *, vector<double> &);

/* IMPLEMENTATIONS */

// Take one time step.
void intg_timeStep(bool sprs, int Lid, int dim, int ord, int var, int n, int M, double h, double k, double *p, double *u0, double *u1){
	int m;
	double u1i[M], u1f[M];

	// initialize
	for (m = 0; m < M; m++){
		u1i[m] = u0[m];
		u1f[m] = u0[m];
	}

	// Newton iteration
	newt_iter(sprs, Lid, dim, ord, var, n, M, h, k, p, u0, u1i, u1f);

	// copy solution to output
	for (m = 0; m < M; m++){
		u1[m] = u1f[m];
	}
}

// Evolve system in time. Assume u is an empty vector upon entry.
// At the end of time evolution, u is populated with all entries
// in time tn and space xm, indexed by u[n*M + m].
void intg_timeEvol(bool sprs, int Lid, int dim, int ord, int var, int N, int M, double h, double k, double *p, double *u0, vector<double> &u){
	int n, m;
	double u1[M];

	for (n = 0; n < N; n++){
		cout << "ts = " << n << " / " << N << endl;
		for (m = 0; m < M; m++)
			u.push_back(u0[m]);
		
		intg_timeStep(sprs, Lid, dim, ord, var, n, M, h, k, p, u0, u1);

		// NEED TO IMPLEMENT REGRIDDING SCHEME FOR THE CASE WHEN THE 
		// TOTAL VOLUME OF A SPREADING FILM IS FIXED (FREE BOUNDARY AT OUTER RIM, x = R(t))
		// IF u0[m] =     gsl_sf_cos(M_PI*x/2.0);
		// THEN Q = 2.0*(M_PI -2.0)/(M_PI*M_PI);
		// REGRIDDING SHOULD BE
//		if (Lid > 1){
//			double sum = 0.0;
//			for (m = 1; m < M-1; m++){
//				sum += m*u1[m];//
//			}
//			dx = sqrt(Q/sum);
////			cout << dx << endl;
//		}

		for (m = 0; m < M; m++){
			u0[m] = u1[m];
		}
	}
}

#endif
