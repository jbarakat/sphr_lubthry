/* MAIN PROGRAM */

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_bessel.h>
#include "../include/diff.h"
#include "../include/syst.h"

using namespace std;

int main(){
	int i, j, k, l;
	int m = 10000;
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

	return 0;
}
