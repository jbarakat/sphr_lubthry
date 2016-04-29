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
#include "../include/syst.h"
#include "../include/diff.h"
#include "../include/newt.h"

using namespace std;

int main(){
	int i, j, k, l;
	int m, n;
	int N = 10;
	int M = 20;
	double p[1];	 // parameters
	double g0[M];	 // source
	double g1[M];	 // source
	double u[M*N]; // function
	bool sprs = true;
	
	double dx = 1.0/double(M-1);
	double dt = dx/6.0;

	double t[N], x[M];
	double u0[M], u1[M];
	double F[M], DF[M*M];

	// space and time grids
	for (m = 0; m < M; m++){
		x [m] = m*dx;
	}

	for (n = 0; n < N; n++)
		t [n] = n*dt;
	
	p[0] = 1.0; // diffusion coefficient

	// Test problem: Heat equation w/Dirichlet BCs.
	int Lid = 0; // equation type
	int o  = 2; // order of spatial operator
	double ui[M], uf[M];

	// initial condition
	for (m = 0; m < M; m++)
		u0[m] =     gsl_sf_sin(M_PI*x[m]);

	// time advance
	for (n = 0; n < N; n++){
		u[n*N + m] = u0[m]; 

		for (m = 0; m < M; m++){
			g0[m] = ((M_PI*M_PI) - 1.0)*gsl_sf_exp(-t[n  ])*gsl_sf_sin(M_PI*x[m]);
			g1[m] = ((M_PI*M_PI) - 1.0)*gsl_sf_exp(-t[n+1])*gsl_sf_sin(M_PI*x[m]);
			ui[m] = u0[m];
			uf[m] = ui[m];
		}

		newt_iter(sprs, Lid, o, n, M, dx, dt, p, g0, g1, u0, ui, uf);

		for(m = 0; m < M; m++){
			u0[m] = uf[m];
		}
	}

	return 0;
}








void testGravitySpreading(){
	int i, j, k, l;
	int m, n;
	int N = 10;
	int M = 100;
	
	/*---- GRAVITATIONAL SPREADING ----*/
	/* SCALINGS:
	 *  h(r,t) scaled by H0 (initial drop height)
	 *  R(t)   scaled by R0 (initial drop radius)
	 *  r      scaled by R0
	 *  t      scaled by mu*R0^2/(rho*g*H0^3) where
	 *    rho = liquid density
	 *    mu  = liquid viscosity
	 *    g   = acceleration due to gravity
	 *  With these scalings, no parameters
	 *  appear in the problem (self-similar).
	 */
	int id = 1; // gravitational spreading
	double r[M], t[N];
	double h[M*N], hi[M];
	double h0[M], h1[M];
	double si[M], sf[M];

	// set parameter to unity (parameter-free problem)
	double p[1];
	p[0] = 1.0;

	// initial drop radius = unity
	double R0 = 1.0;
	
	// time and grid step sizes
	// (grid step size will change at each timestep
	//  in order to conserve volume)
	double dr = R0/double(M-1);
	double dt = dr/6.0;

	// time domain discretization
	for (n = 0; n < N; n++){
		t[n] = n*dt;
	}

	// space domain discretization (radial coordinate)
	for (m = 0; m < M; m++){
		r[m] = m*dr;
	}
	
	// initial thickness profile
	for (m = 0; m < M; m++){
		hi[m] = gsl_sf_cos(M_PI*r[m]/2.0);
	}

	// TEST NEWTON SOLVER USING INITIAL PROFILE
	for (m = 0; m < M; m++){
		h0[m] = hi[m];
		h1[m] = hi[m];
		
		// guess s for h1
		si[m] = h0[m];
		sf[m] = h0[m];
	}

	//newt_iter(id, M, dr, dt, p, h0, si, sf);

	// update dr
	// (will have to have clever algorithm here
	//  to ensure that the grid points are not
	//  too widely spaced apart)
	
}








void testSyst(){
	int i, j, k, l;
	int m, n;
	int N = 10;
	int M = 10;
	
	double dr = 2*M_PI/double(M-1);
	double dt = dr/6.0;

	double t[N], r[M];
	double u0[M], u1[M];
	double F[M], DF[M*M];
	for (i = 0; i < M; i++){
		r [i] = i*dr;
		u0[i] =      gsl_sf_cos(r[i]);
		u1[i] = 1.1*gsl_sf_cos(r[i]);
	}
	
	int id = 0; // equation type
	int o  = 2; // order of spatial operator
	double p[1];// parameters
	double g[M];// source
	bool sprs = true;
	p[0] = 1.0;
	for (m = 0; m < M; m++) g[m] = 0.0;
//	syst_eqn(id, M, dr, dt, p, u0, u1, F, DF);
	syst_F_DF(sprs, id, o, n, M, dr, dt, p, g, g, u0, u1, F, DF);

	for (i = 0; i < M; i++){
		for (j = 0; j < M; j++){
			if (DF[i*M + j] < 0)
				printf( "%.4f ", DF[i*M+j]);
			else
				printf(" %.4f ", DF[i*M+j]);
		}
		printf("\n");
	}
	
	double si[M], sf[M];
	for (m = 0; m < M; m++){
		si[m] = u1[m];
		sf[m] = 0.0;
	}
	newt_iter(sprs, id, o, n, M, dr, dt, p, g, g, u0, si, sf);
}

void testDiff(){
	int i, j, k, l;
	int m = 1000;
	vector<double> t(m), cost(m), sint(m);
	vector<double> dsint(m), dcost(m);

	vector<double> x(m), K0(m), I0(m);
	vector<double> bK0(m), bI0(m);

	for (i = 0; i < m; i++){
		x[i] = i*10.0/m;
		t[i] = i*M_PI/m;
		cost[i] = gsl_sf_cos(t[i]);
		sint[i] = gsl_sf_sin(t[i]);
		//K0  [i] = gsl_sf_bessel_K0(x[i]);
		I0  [i] = gsl_sf_bessel_I0(x[i]);
	}

	double dt = t[1] - t[0];
	double dx = x[1] - x[0];

	diff_d4(m, dt, sint.data(), dsint.data());
	diff_b2(m, dx, I0.data(), bI0.data());

	double error = 0.0;
	for (i = 0; i < m; i++){
		//error += fabs(dsint[i] + cost[i]);
		//error += fabs(dsint[i] - sint[i]);

		//cout << fabs(dsint[i] + cost[i]) << endl;
		//cout << fabs(dsint[i] - sint[i]) << endl;

		error += fabs(bI0[i]);
		if (i < 10)
			cout << bI0[i] << endl;
	}
	error /= m;

	cout << error << endl;

}
