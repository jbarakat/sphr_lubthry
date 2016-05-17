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

void init(int, int, int, double, double *, double *);

int main(){
	int i, j, k, l;
	int m, n;
	
	// parameters
	double p[10];						 	// parameters
				 p[0] = 1.0;			 	// diffusion coefficient
				 p[1] = 1.0;				// concavity coefficient (for lub thry on curved surface)
	int sprs = 0;					// sparsity indicator
	string dir = "../output";	// output directory
	string fn = "data";				// output file name
	int nw = 4 ;							// number of timesteps to write

	int dim = 1;						 	// dimensions in space
	int var = 2;							// number of dependent variables
	
	// setup time and space domains
	int N = 40;										// number of time points
	int M = 100;									// number of grid points (even)
	int Mv = var*M;								// size of system
	double xmax = 2.0;
	double tmax = 0.01;
	double dx = xmax/double(M-1);
	double dt = tmax/double(N-1);



	// NUMBER OF VARIABLES = 1 (diffusion of one species)

	// NUMBER OF VARIABLES = 2
	// Make film thickness and flux separate variables
	// eventually have surface concentration be a separate variable too...
	// for now... just double # of grid points...

	double u0[Mv];						// initial condition
	vector<double> u;					// solution vector

	// initialize
	init(dim, var, M, dx, p, u0);

	// evolve system in time
	intg_timeEvol(sprs, dim, var, N, M, dx, dt, p, u0, u);

	// write to file
	write(dir, fn, nw, var, N, M, dx, dt, u.data());

	return 0;
}

void init(int dim, int var, int M, double h, double *p, double *u0){
	int m, v;
	double x;
	double cf1 = p[0];
	double H[M], Q[M];
	double H2, H3, d3H;
	
	// initialize H
	for (m = 0; m < M; m++){
		x = m*h;
		if (x < 1.0) H[m] = 1.0 - 4.0*(x - 0.5)*(x - 0.5);
		else         H[m] = 0.1;
	}

	for (m = 0; m < M; m++){
		x = m*h;
		H2 = H[m]*H[m];
		H3 = H2*H[m];
		diff_d3(m, M, h, H, d3H);
		Q[m]  = (1.0/3.0)*H3*(1.0 + cf1*d3H);
		u0[var*m + 0] = H[m];
		u0[var*m + 1] = Q[m];
	}
}

