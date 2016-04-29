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
void intg_timeStep(bool sprs, int Lid, int gid, int o, int n, int M, double h, double k,
                   double *p, double *u0, double *u1){
	int m;
	double u1i[M], u1f[M];
	double g0 [M], g1 [M];
	double g0m   , g1m   ;

	// initialize
	for (m = 0; m < M; m++){
		u1i[m] = u0[m];
		u1f[m] = u0[m];
		
		prbm_g(gid, n  , m, M, h, k, p, g0m);
		prbm_g(gid, n+1, m, M, h, k, p, g1m);
		
		g0[m] = g0m;
		g1[m] = g1m;
	}

	// Newton iteration
	newt_iter(sprs, Lid, o, n, M, h, k, p, g0, g1, u0, u1i, u1f);

	// copy solution to outptu
	for (m = 0; m < M; m++){
		u1[m] = u1f[m];
	}
}

#endif
