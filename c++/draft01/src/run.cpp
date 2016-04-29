/* MAIN PROGRAM */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include <chplot.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_bessel.h>
#include "../include/intg.h"

using namespace std;

int main(){
	int i, j, k, l;
	int m, n;
	int N = 100;
	int M = 20;
	double p[1];	 // parameters
				 p[0] = 1.0; // diffusion coefficient
	vector<double> u;
	bool sprs = true;
	double dx = 1.0/double(M-1);
	double dt = dx/6.0;

	double u0[M], u1[M];
	

	// Test problem: Heat equation w/Dirichlet BCs.
	int Lid = 0; // operator indicator
	int gid = 1; // source indicator
	int o  = 2; // order of spatial operator

	// initial condition
	for (m = 0; m < M; m++)
		u0[m] =     gsl_sf_sin(M_PI*m*dx);

	// time advance
	for (n = 0; n < N; n++){
		u.push_back(u0[m]);

		intg_timeStep(sprs, Lid, gid, o, n, M, dx, dt, p, u0, u1);

		for(m = 0; m < M; m++){
			u0[m] = u1[m];
		}
	}

	return 0;
}
