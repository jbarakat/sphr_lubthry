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

using namespace std;

int main(){
	int i, j, k, l;
	int m = 50;
	vector<double> t(m), cost(m), sint(m);
	vector<double> dsint(m), dcost(m);

	vector<double> x(m), K0(m), I0(m);
	vector<double> bK0(m), bI0(m);

	for (i = 0; i < m; i++){
		x[i] = 0.001 + i*10.0/m;
		t[i] = i*M_PI/m;
		cost[i] = gsl_sf_cos(t[i]);
		sint[i] = gsl_sf_sin(t[i]);
		K0  [i] = gsl_sf_bessel_K0(x[i]);
		I0  [i] = gsl_sf_bessel_I0(x[i]);
	}

	double h = t[1] - t[0];

	d4(m, h, sint.data(), dsint.data());
	b2(m, h, 0.001, I0.data(), bI0.data());

	double error = 0.0;
	for (i = 0; i < m; i++){
		//error += fabs(dsint[i] + cost[i]);
		//error += fabs(dsint[i] - sint[i]);

		//cout << fabs(dsint[i] + cost[i]) << endl;
		//cout << fabs(dsint[i] - sint[i]) << endl;

		error += fabs(bI0[i]);
		cout << bI0[i] << endl;
	}
	error /= m;

	cout << error << endl;

	return 0;
}
