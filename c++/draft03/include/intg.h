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
void intg_timeStep(int, int, int, int, int, double, double, double *, double *, double *);
void intg_timeEvol(int, int, int, int, int, double, double, double *, double *, vector<double> &);

/* IMPLEMENTATIONS */

// Take one time step.
void intg_timeStep(int sprs, int dim, int var, int n, int M, double h, double k, double *p, double *u0, double *u1){
	int m;
	int Mv = var*M;
	double u1i[Mv], u1f[Mv];

	// initialize
	for (m = 0; m < Mv; m++){
		u1i[m] = u0[m];
		u1f[m] = u0[m];
	}

	// Newton iteration
	newt_iter(sprs, dim, var, n, M, h, k, p, u0, u1i, u1f);

	// copy solution to output
	for (m = 0; m < Mv; m++){
		u1[m] = u1f[m];
	}
}

// Evolve system in time. Assume u is an empty vector upon entry.
// At the end of time evolution, u is populated with all entries
// in time tn and space xm, indexed by u[n*M + m].
void intg_timeEvol(int sprs, int dim, int var, int N, int M, double h, double k, double *p, double *u0, vector<double> &u){
	int n, m;
	int Mv = var*M;
	double u1[Mv];

	for (n = 0; n < N; n++){
		cout << "ts = " << n << " / " << N << endl;
		for (m = 0; m < Mv; m++)
			u.push_back(u0[m]);
		
		intg_timeStep(sprs, dim, var, n, M, h, k, p, u0, u1);

		for (m = 0; m < Mv; m++){
			u0[m] = u1[m];
		}
	}
}

#endif
