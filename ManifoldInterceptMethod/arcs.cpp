#include "arcs.h"
#include <stdexcept>
#include <limits>

arc::arc() : order(-1), sign(true), tangent(-1), time(-1) {}

arc::arc(int order_, bool sign_, int tangent_, double time_/*=0*/) : order(order_), sign(sign_), tangent(tangent_), time(time_) {}

void Constraint::clear() {
	if (M_max) {
		delete[] M_max;
		M_max = nullptr;
	}
	if (M_min) {
		delete[] M_min;
		M_min = nullptr;
	}
}

// Copy the information of the constraint.
void Constraint::copy(int order, const double* M_max, const double* M_min) {
	if (isinf(M_max[0]) || isinf(M_min[0])) {
		throw invalid_argument("M_max[0] and M_min[0] must be finite.");
		exit(1);
	}
	clear();
	this->M_max = new double[order + 1];
	this->M_min = new double[order + 1];
	for (int i = 0; i <= order; i++) {
		this->M_max[i] = M_max[i];
		this->M_min[i] = M_min[i];
	}
}

Constraint::~Constraint() {
	clear();
}

// Change the sign of the constraint.
void Constraint::change_sign(int order) {
	swap(M_max, M_min);
	for (int i = 0; i <= order; i++) {
		M_max[i] = -M_max[i];
		M_min[i] = -M_min[i];
	}
}

void waypoints::clear() {
	if (xs != nullptr) {
		delete[] xs;
		xs = nullptr;
		order = 0;
		length = 0;
	}
}

void waypoints::set(int order, int length) {
	if (xs == nullptr) {
		xs = new double[order * length]();
	}
	else if ((this->order * this->length) < (order * length)) {
		delete[] xs;
		xs = new double[order * length]();
	}
	this->order = order;
	this->length = length;
}

waypoints::~waypoints() {
	clear();
}

waypoints::waypoints() : order(0), length(0), xs(nullptr) {}