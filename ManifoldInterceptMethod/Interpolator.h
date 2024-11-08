#pragma once
#include <vector>
#include <fstream>
#include "arcs.h"
using namespace std;

class Interpolator {
public:
	Interpolator(int order);
	~Interpolator();
	double interpolate(const double* x0, const arc* arcs, int length, const Constraint& constraint, double Ts, double T0 = 0, bool flag_end = true);
	void write_csv(const string& filename);
public:
	vector<double>* buffer;
	int order;
};