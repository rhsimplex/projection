// projection.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include "armadillo"
#include "E8projection.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
	string inputfile = argv[1];
	int numpoints;
	cout << "Started: " << __DATE__ << " " << __TIME__ << endl;
	E8projection Q(inputfile);
	Q.printprojectionmatrix("Projection Matrix: ");
	Q.printbasis3D("3D Basis: ");
	Q.printshiftvector("Shift Vector: ");
	numpoints = Q.execute(inputfile.append(".out"));
	cout << "Number of points projected: " << numpoints << endl;
}

