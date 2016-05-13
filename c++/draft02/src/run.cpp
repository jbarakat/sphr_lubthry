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

void init(int, int, int, int, int, double, double *);

int main(){
	int i, j, k, l;
	int m, n;
	
	// parameters
	double p[10];						 	// parameters
				 p[0] = 1.0;			 	// diffusion coefficient
				 p[1] = 1.0;				// concavity coefficient (for lub thry on curved surface)
	bool sprs = true;					// sparsity indicator
	string dir = "../output";	// output directory
	string fn = "data";				// output file name
	int nw = 4 ;							// number of timesteps to write

	// CHANGE THE PROBLEM TYPE... NEED TO BE MORE CAREFUL
	// WHEN I GET TO THE REAL PROBLEM...
	int Lid = 1;							// refer to prbm.h for ID definitions
	int dim = 1;						 	// dimensions in space
	int var = 2;							// number of dependent variables
	int ord = 4;							// order of spatial operator
	
	// setup time and space domains
	int N = 4;										// number of time points
	int M = 8;										// number of grid points (even)
	double xmax = 4.0;
	double tmax = 0.01;
	double dx = xmax/double(M/var-1);
	double dt = tmax/double(N-1);


	// NUMBER OF VARIABLES = 1 (diffusion of one species)

	// NUMBER OF VARIABLES = 2
	// Make film thickness and flux separate variables
	// eventually have surface concentration be a separate variable too...
	// for now... just double # of grid points...

	double u0[M];							// initial condition
	vector<double> u;					// solution vector

	// initialize
	init(Lid, dim, ord, var, M, dx, u0);

	// evolve system in time
	intg_timeEvol(sprs, Lid, dim, ord, var, N, M, dx, dt, p, u0, u);

	// write to file
	write(dir, fn, nw, var, N, M, dx, dt, u.data());

	return 0;
}

void init(int Lid, int dim, int ord, int var, int M, double h, double *u0){
	int m;
	double x;
	double mu, sig, sig2;
	sig  = 0.05;
	sig2 = sig*sig;

	for (m = 0; m < M; m++){
		x = m*h;
		u0[m] = 0.0;
		if (var == 1){
			if (ord == 2){
				if (Lid == 0){	// diffusion of Gaussian source
					if (dim == 1) mu = 0.5;
					if (dim == 2) mu = 0.0;
					u0[m] = gsl_sf_exp(-(x - mu)*(x - mu)/sig2)/(sqrt(2.0*M_PI*sig2));
					cout << u0[m] << endl;
				}
				if (Lid == 1){	// gravitational spreading of drop
					if (dim == 1) u0[m] = gsl_sf_sin(M_PI*x    );
					if (dim == 2) u0[m] = gsl_sf_cos(M_PI*x/2.0);
				}
			}
			if (ord == 4){
				if (Lid == 2){	// penetration of rigid sphere through flat interface
					u0[m] = 0.0;
				}
			}
		}

		if (var == 2){
			if (ord == 4){
				int m2;
				int M2 = M/2;
				double H, Q;
				if (Lid == 1){	// gravity/capillary spreading
					if (dim == 1){
						if      (m < M2){	// film thickness H(x,t)
							m2 = m;
							x  = m2*h;
							if (x < 1) H  = 1.0 - 4.0*(x - 0.5)*(x - 0.5);
							else       H  = 0.1;
							u0[m] = H;
						}
						else            { // flux Q(x,t)
							m2 = m - M2;
							x  = m2*h;
							if (x < 1) H  = 1.0 - 4.0*(x - 0.5)*(x - 0.5);
							else       H  = 0.1;
							double H2, H3, d3H;
							H2 = H*H;
							H3 = H2*H;
							diff_d3(m-M2, M2, h, u0, d3H);
							double PAR = 1.0; // NEED TO EDIT THIS... PASS PARAMETER TO INITIAL CONDITION FUNCTION...
							Q  = (1.0/3.0)*H3*(1.0 + PAR*d3H);
							u0[m] = Q;
						}
					}
				}
			}
		}
	}
}

