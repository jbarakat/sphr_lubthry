/* WRITE FILE
 *  Write data file containing the abscissa and solution vectors.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 */

#ifndef WRITE_H
#define WRITE_H

/* HEADER FILES */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

/* PROTOTYPES */
void write(string, string, int, int, int, double, double, double *);
void write(string, string, int, int, int, int, double, double, double *);

/* IMPLEMENTATIONS */

// Write output to file.
void write(string dir, string fn, int var, int n, int M, double h, double k, double *uw){
	int m;
	int width = 12;
	int tacc = 4;
	int xacc = 4;
	int uacc = 4;
	int fill = 6;
	string fmt = ".txt";
	string col1 = "x";
	string col2 = "u";
	ofstream file;
	ostringstream time, ts, params, header;
	string line, path;

	// setup: assign strings and stringstreams
	time << fixed << setprecision(tacc) << n*k;
	ts << setfill('0') << setw(fill) << n;
	params << "t = " << time.str();
	cout << params.str() << endl;
	header << setw(width) << col1 << setw(width) << col2;
	path = dir + "/" + fn + ts.str() + fmt;

	// write to file
	file.open(path.c_str());
	file << params.str() << endl;
	file << header.str() << endl;
	if (var == 1){
		for (m = 0; m < M; m++){
			ostringstream xx, uu;
			xx << fixed << setprecision(xacc) << setw(width) << m*h;
			uu << fixed << setprecision(uacc) << setw(width) << uw[m];
			
			line = xx.str() + uu.str();
			file << line << endl;
		}
	}
	else if (var == 2){
		int M2 = M/var;	
		for (m = 0; m < M2; m++){
			ostringstream xx, u1, u2;
			xx << fixed << setprecision(xacc) << setw(width) << m*h;
			u1 << fixed << setprecision(uacc) << setw(width) << uw[m];
			u2 << fixed << setprecision(uacc) << setw(width) << uw[m+M2];
			
			line = xx.str() + u1.str() + u2.str();
			file << line << endl;
		}
	}
	else if (var == 3){
		int M3 = M/var;	
		for (m = 0; m < M3; m++){
			ostringstream xx, u1, u2, u3;
			xx << fixed << setprecision(xacc) << setw(width) << m*h;
			u1 << fixed << setprecision(uacc) << setw(width) << uw[m];
			u2 << fixed << setprecision(uacc) << setw(width) << uw[m+M3];
			u3 << fixed << setprecision(uacc) << setw(width) << uw[m+2*M3];
			
			line = xx.str() + u1.str() + u2.str() + u3.str();
			file << line << endl;
		}
	}
	file.close();

	cout << "Output written to " << path << "." << endl;
}

void write(string dir, string fn, int nw, int var, int N, int M, double h, double k, double *u){
	int n, m;
	double uw[M];

	for (n = 0; n < N; n++){
		if (n % nw == 0){
			for (m = 0; m < M; m++)
				uw[m] = u[n*M + m];
			write(dir, fn, var, n, M, h, k, uw);
		}
		else
			continue;
	}
}


#endif
