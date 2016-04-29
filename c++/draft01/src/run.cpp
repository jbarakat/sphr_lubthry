/* MAIN PROGRAM */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_bessel.h>
#include "../include/intg.h"
#include "../include/write.h"

using namespace std;

int main(){
	int i, j, k, l;
	int m, n;
	
	/*----------- SETUP -----------*/
	
	// setup time and space domains
	int N = 40;
	int M = 100;
	double xmax = 1.0;
	double tmax = 0.01;
	double dx = xmax/double(M-1);
	double dt = tmax/double(N-1);

	// parameters
	double p[1];	 // parameters
				 p[0] = 1.0; // diffusion coefficient
	
	// exploit sparsity in iterative solve?
	bool sprs = true;

	// solution vector	
	vector<double> u;
	double u0[M], u1[M];
	
	/*----------- SOLVE -----------*/

	// Test problem: Heat equation w/Dirichlet BCs.
	int Lid = 1; // operator type
	int o   = 2; // order of spatial operator
	double mu  = 0.5;
	if (Lid > 0)
		mu = 0.0;
	double sig = 0.05;

	// initial condition
	for (m = 0; m < M; m++){
		double x = m*dx;

		// diffusion problems: intial profile is sine wave or Gaussian
		if (Lid == 0 || Lid == 1){
			u0[m] =     gsl_sf_sin(M_PI*x);
			u0[m] =     gsl_sf_exp(-(x - mu)*(x - mu)/(sig*sig))/(sqrt(2.0*M_PI*sig*sig));
		}

		// spreading problems: initial profile is a cosine
	}

	// time advance
	for (n = 0; n < N; n++){
		for (m = 0; m < M; m++){
			u.push_back(u0[m]);
		}

		intg_timeStep(sprs, Lid, o, n, M, dx, dt, p, u0, u1);

		for(m = 0; m < M; m++){
			u0[m] = u1[m];
		}
	}
	
	/*-------- WRITE TO FILE -------*/
	string dir = "../output";
	string fn = "data";

	for (n = 0; n < N; n++){
		if (n % 10 == 0){
			for (m = 0; m < M; m++){
				u0[m] = u[n*M + m];
			}
			write(dir, fn, n, M, dx, dt, u0);
		}
		else
			continue;
	}

	return 0;
}
