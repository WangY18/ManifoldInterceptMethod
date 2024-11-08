#include "Interpolator.h"
#include "Planner.h"
#include <iostream>

Interpolator::Interpolator(int order_): order(order_) {
	buffer = new vector<double>[order + 1]();
}

Interpolator::~Interpolator() {
	if (buffer != nullptr) {
		delete[] buffer;
	}
}

// Interpolate the states and control. The initial time is T0. 
// T0: the initial time. Ts: the sample time. flag_end: whether the last arc is the end point.
// Return: the remaining time. (If flag_end is true, the remaining time is 0.)
double Interpolator::interpolate(const double* x0, const arc* arcs, int length, const Constraint& constraint, double Ts, double T0/*=0*/, bool flag_end/*=true*/) {
	double* x0now = new double[order]();
	double* xnow = new double[order]();
	for (int k = 0; k < order; k++) {
		x0now[k] = x0[k];
	}

	bool flag_interpoalte = false;
	double u = 0;
	for (int i = 0; i < length; i++) {
		if (arcs[i].tangent != 0) {
			continue;
		}
		u = constraint.cal_u(arcs[i]);
		while (T0 <= arcs[i].time) {
			Planner::dynamics_onestep(order, x0now, xnow, u, T0);
			buffer[0].push_back(u);
			for (int k = 0; k < order; k++) {
				buffer[k + 1].push_back(xnow[k]);
			}
			T0 += Ts;
			flag_interpoalte = true;
		}
		T0 -= arcs[i].time;
		Planner::dynamics_onestep(order, x0now, xnow, u, arcs[i].time);
		swap(xnow, x0now);
	}

	if (flag_end && fabs(T0 - Ts) > 1e-6) {
		buffer[0].push_back(u);
		for (int k = 0; k < order; k++) {
			buffer[k + 1].push_back(x0now[k]);
		}
	}

	delete[] x0now;
	delete[] xnow;
	return Ts - T0;
}

// Output the states and control to a csv file.
void Interpolator::write_csv(const string& filename) {
	std::ofstream file(filename);
	if (!file.is_open()) {
		cout << "Error opening file: " << filename << std::endl;
		return;
	}
	file << "u";
	for (int k = 1; k <= order; ++k) {
		file << ",x" << k;
	}
	file << "\n";
	for (int i = 0; i < buffer[0].size(); ++i) {
		file << buffer[0][i];
		for (int k = 1; k <= order; ++k) {
			file << "," << buffer[k][i];
		}
		file << "\n";
	}
	file.close();
}