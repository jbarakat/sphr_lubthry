/* NEWTON-RAPHSON METHOD
 *  Solve the linear system DF*d = F for the Newton step d.
 *
 * REFERENCES
 *  Stoer & Bulirsch, pp. 302-316, 559-561
 *  
 * PARAMETERS
 *  n		[input]			number of ODEs
 *  m		[input]			number of shooting points
 *  t		[input]			independent coordinate (ranges from 0 to 1)
 *  s		[input]			guess (an mn-vector)
 *  y		[input]			solution [an (m-1)n-vector]
 *  F		[output]		function to be zeroed
 *  DF	[output]		Jacobian
 */

#ifndef NEWTON_H
#define NEWTON_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include "./syst.h"

using namespace std;

/* PROTOTYPES */
void newt_d(int, int, double, double, double *, double *, double *, double *, double *, double*);

/* IMPLEMENTATIONS */

// Iterate on f1 until F(f1;f0,p) = 0
void newt_iter(int ID, int M, double h, double k, double *p, double *f0, double *f1){
	int    i, j, m;
	double f1k[M];
	double F[M], DF[M], d[M];
	double Fk[M], DFk[M], dk[M];

	const int MAXITER = 1000;
	const double TOL  = 1e-6;

	// initialize
	int    iter  = 0;			// iteration counter
	double Fnorm = 1.0;		// 2-norm of F
	double nu    = 1.0;		// fraction of full Newton step

	// iterate until F(f1;f0,p) = 0
	while (iter < MAXITER && fnorm > TOL){
		// generate F, DF, and d
		newt_d(ID, M, h, k, p, f0, f1, F, DF, d);

		// calculate 2-norm of F
		Fnorm = 0;
		for (m = 0; m < M; m++)
			Fnorm += F[m]*F[m];
		Fnorm = sqrt(Fnorm);

		// use a finite search process to determine
		// the fraction of the full Newton step
		for (i = 0; i < 5; i++){
			nu = pow(2.0, -i);

			for (m = 0; m < M; m++){
				f1k[m] = f1[m] - nu*d[j];
			}
			
			// recalculate Fk using f1k
			newt_d(ID, M, h, k, p, f0, f1k, Fk, DFk, dk);
			
			Fknorm = 0;
			for (m = 0; m < M; m++)
				Fknorm += Fk[m]*Fk[m];
			Fknorm = sqrt(Fknorm);
			
			if (Fknorm < Fnorm)
				break;
		}

		// update f1
		for (m = 0; m < M; m++)
			f1[m] = f1k[m];
	}
}


// Compute the full Newton step d, where DF*d = F
//  f0 = solution at previous timestep
//  f1 = guess for solution at new timestep
void newt_d(int ID, int M, double h, double k, double *p, double *f0, double *f1, double *F, double *DF, double *d){
	int i, j, info;
	double Fp[M], DFp[M*M];
	
	// given f0 and f1, generate F (function) and DF (Jacobian)
	syst_eqn(ID, M, h, k, p, f0, f1, Fp, DFp);

	// copy F and DF to output
	// (Fp will be transformed by LAPACK operations)
	for (i = 0; i < M; i++){
		F[i] = Fp[i];
		for (j = 0; j < M; j++){
			DF[i*M + j] = DFp[i*M + j];
		}
	}
	
	/* solve the linear system using Gaussian elimination 
	 * (LAPACK's DGETRF and DGETRS subroutines) */
	char   trans = 'N';				// form of the system of equations
	int    ldDF  =  M;				// leading dimension of DF
	int    ldF   =  1;				// leading dimension of F
	int    nrhs  =  1;				// number of right-hand sides
	int    ipiv[M];		 				// array of pivot indices

	// factorize
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, M, M, DFp, ldDF, ipiv);
	if (info != 0){
		cout << "Error: dgetrf did not return successfully." << endl;
		return;
	}

	// solve
	info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, trans, M, nrhs, DFp, ldDF, ipiv, Fp, ldF);
	if (info != 0){
		cout << "Error: dgetrs did not return successfully." << endl;
		return;
	}
	
	// store solution
	for (i = 0; i < M; i++){
		d[i] = Fp[i];
	}

//	//// LATER, IMPLEMENT BROYDEN'S UPDATE
//	//// -- note that Broyden shouldn't be used too many times in sequence
//	//// -- near the solution, Broyden's method becomes inaccurate
//	//// -- forward differences miiiight be inaccurate near the solution as well..
//	  //    (will need to use central difference)
}

#endif
