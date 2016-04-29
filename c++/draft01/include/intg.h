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
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include "./newt.h"

/* PROTOTYPES */

/* IMPLEMENTATIONS */
void intg_timeStep(bool sprs, int Lid, int o, int n, int M, double h, double k,
                   double *p, double *u0, double *u1){
	int m;
	double u1i[M], u1f[M];

	// initialize
	for (m = 0; m < M; m++){
		u1i[m] = u0[m];
		u1f[m] = u0[m];
	}

	// Newton iteration
	newt_iter(sprs, Lid, o, n, M, h, k, p, u0, u1i, u1f);

	// copy solution to outptu
	for (m = 0; m < M; m++){
		u1[m] = u1f[m];
	}
}

#endif
