#pragma once
#include "armadillo"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace arma;
using namespace std;
class E8projection
{
public:
	E8projection(string filename);
	E8projection(float t);
	void project(colvec point8D);
	~E8projection(void);
	void printprojectionmatrix(string name);
	void printprojectionmatrix();
	void printbasis3D(string name);
	void printbasis3D();
	void printorthonormalbasis8D(string name);
	void printorthonormalbasis8D();
	void printshiftvector(string name);
	int execute(string output);
	static mat orthonormalize(const mat &toorth);
	colvec current3Dpoint;
	colvec current3Dpointcc;
	colvec current8Dpoint;
	float currentprojdist;
protected:
	mat E8D3D;
	mat E8D3Dtr;
	mat orthonormalbasis3D;
	mat basis3D;
	mat scale3D;
	colvec shiftvector;
	int LAYERS;
	float MAXdist;
private:
	ifstream inputfile;
};

