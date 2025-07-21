#include "Planner_MIM.h"
#include <gsl/gsl_poly.h>
#include <stdexcept>
#include <iostream>


namespace MIM {

#include <cmath>

// Plan a 1st-order trajectory_SCP from x0 to xf under constraint. (trivial)
inline arc Planner_MIM::plan_1st_order(double x0, double xf, const Constraint& constraint) {
	return x0 < xf ? arc(0, true, 0, (xf - x0) / constraint.M_max[0]) : arc(0, false, 0, (xf - x0) / constraint.M_min[0]);
}

// Plan a 2nd-order trajectory_SCP from x0 to xf under constraint.
// Only: 00, 010
// If direction is 1/-1, then directly consider the case where x0 is lower/higher than xf.
vector<arc> Planner_MIM::plan_2nd_order(double* x0, double* xf, Constraint& constraint, bool flag_consider_position, int direction) {
	if (direction == 0) {
		arc arc_try = plan_1st_order(x0[0], xf[0], constraint);
		double* xs_try = end_point(2, x0, constraint, &arc_try, 1);
		if (norm_inf(xs_try, xf, 2) < constraint.EPSILON) {
			vector<arc> result;
			result.push_back(arc_try);
			if (flag_consider_position && (feasible(2, x0, constraint, result.data(), result.size()) != 0)) {
				result.clear();
			}
			delete[] xs_try;
			return result;
		}
		direction = xs_try[1] < xf[1] ? 1 : -1;
		delete[] xs_try;
	}
	if (direction < 0) { // x0 is higher than xf. It must be feasible.
		constraint.change_sign(2);
		change_sign_x(x0, 2);
		change_sign_x(xf, 2);
		vector<arc> result = plan_2nd_order(x0, xf, constraint, flag_consider_position, 1);
		constraint.change_sign(2);
		change_sign_x(x0, 2);
		change_sign_x(xf, 2);
		change_sign_arcs(result.data(), result.size());
		return result;
	}
	// x0 is lower than xf
	vector<arc> result;
	double DeltaX2 = xf[1] - x0[1] + (constraint.M_max[1] * constraint.M_max[1] - xf[0] * xf[0]) / (2.0 * constraint.M_min[0]) - (constraint.M_max[1] * constraint.M_max[1] - x0[0] * x0[0]) / (2.0 * constraint.M_max[0]);
	if ((!isinf(constraint.M_max[1])) && DeltaX2 > 0) { // Consider 010
		result.push_back(arc(0, true, 0, (constraint.M_max[1] - x0[0]) / constraint.M_max[0]));
		result.push_back(arc(1, true, 0, DeltaX2 / constraint.M_max[1]));
		result.push_back(arc(0, false, 0, (xf[0] - constraint.M_max[1]) / constraint.M_min[0]));
		if (flag_consider_position && feasible(2, x0, constraint, result.data(), result.size()) != 0) {
			result.clear();
		}
		return result;
	}
	// Consider 00
	// coeffs[0] + coeffs[1]* T2 + coeffs[2] * T2^2 = 0
	double coeffs_00[3] = { xf[0] * xf[0] - x0[0] * x0[0] + 2.0 * constraint.M_max[0] * (x0[1] - xf[1]),2.0 * xf[0] * (constraint.M_max[0] - constraint.M_min[0]),constraint.M_min[0] * (constraint.M_min[0] - constraint.M_max[0]) };
	static double x_00[2];
	int num_root = solve_polynomial_real(2, coeffs_00, x_00);
	double T2 = num_root > 1 ? max(x_00[0], x_00[1]) : x_00[0]; // T2 > 0 holds

	double T1 = (xf[0] - x0[0] - (constraint.M_min[0] - constraint.M_max[0]) * T2) / constraint.M_max[0]; // T1 > T2 holds
	result.push_back(arc(0, true, 0, T1 - T2));
	result.push_back(arc(0, false, 0, T2));
	if (flag_consider_position && feasible(2, x0, constraint, result.data(), result.size()) != 0) {
		result.clear();
	}
	return result;
}

// Change the sign of x.
inline void Planner_MIM::change_sign_x(double* x, int order) {
	for (int k = 0; k < order; k++) {
		x[k] = -x[k];
	}
}

inline void Planner_MIM::change_sign_arcs(arc* arcs, int length) {
	for (int i = 0; i < length; ++i) {
		arcs[i].sign = !arcs[i].sign;
	}
}

vector<arc> Planner_MIM::plan(int order, double* x0, double* xf, Constraint& constraint, bool flag_consider_position) {
	if (order == 1) {
		vector<arc> result;
		result.push_back(plan_1st_order(x0[0], xf[0], constraint));
		return result;
	}
	if (order == 2) {
		return plan_2nd_order(x0, xf, constraint, flag_consider_position, 0);
	}
	if (order == 3) {
		return plan_3rd_order(x0, xf, constraint, flag_consider_position);
	}
	return vector<arc>(); // TODO: higher order
}

// Calculate the state vector at the end of the arc.
void Planner_MIM::dynamics_onestep(int order, const double* x0, double* xf, double u, double dt) {
	double T = 1;
	for (int k = 0; k < order; k++) {
		xf[k] = 0;
	}
	for (int i = 0; i < order; i++) {
		// T = dt^k/k!
		for (int k = i; k < order; k++) {
			xf[k] += x0[k - i] * T;
		}
		T *= dt / (i + 1);
		xf[i] += u * T;
	}
}

// Calculate the end points of the trajectory_SCP from the initial state vector x0.
waypoints Planner_MIM::end_points(int order, const double* x0, const Constraint& constraint, const arc* arcs, int length) {
	int N_real = 0;
	for (int i = 0; i < length; i++) {
		if (arcs[i].tangent == 0) {
			N_real++;
		}
	}
	// result[i*order+k] is x[i][k].
	waypoints result;
	result.set(order, N_real + 1);
	for (int k = 0; k < order; k++) {
		result[0][k] = x0[k];
	}
	for (int i = 0, i_arc = 0; i < N_real; i++, i_arc++) {
		while (arcs[i_arc].tangent != 0) {
			i_arc++;
		}
		// x moves from x[k] to x[k+1] along arcs[i_arc]
		dynamics_onestep(order, result[i], result[i + 1], constraint.cal_u(arcs[i_arc]), arcs[i_arc].time);
	}
	return result;
}

// Calculate the end points of the trajectory_SCP from the initial state vector x0.
double* Planner_MIM::end_point(int order, const double* x0, const Constraint& constraint, const arc* arcs, int length) {
	double* x1 = new double[order];
	double* xf = new double[order];
	for (int k = 0; k < order; k++) {
		xf[k] = x0[k];
	}
	for (int i = 0; i < length; i++) {
		if (arcs[i].tangent == 0) {
			dynamics_onestep(order, xf, x1, constraint.cal_u(arcs[i]), arcs[i].time);
			swap(x1, xf);
		}
	}
	delete[] x1;
	return xf;
}

// Calculate the start points of the trajectory_SCP to the terminal state vector xf.
waypoints Planner_MIM::start_points(int order, const double* xf, const Constraint& constraint, const arc* arcs, int length) {
	arc* arcs_reverse = new arc[length];
	for (int i = 0; i < length; i++) {
		arcs_reverse[i] = arcs[length - i - 1];
		arcs_reverse[i].time = -arcs_reverse[i].time;
	}
	waypoints result = end_points(order, xf, constraint, arcs_reverse, length);
	delete[] arcs_reverse;
	return result;
}

// Calculate the end points of the trajectory_SCP from the initial state vector x0.
double* Planner_MIM::start_point(int order, const double* xf, const Constraint& constraint, const arc* arcs, int length) {
	arc* arcs_reverse = new arc[length];
	for (int i = 0; i < length; i++) {
		arcs_reverse[i] = arcs[length - i - 1];
		arcs_reverse[i].time = -arcs_reverse[i].time;
	}
	double* x0 = end_point(order, xf, constraint, arcs_reverse, length);
	delete[] arcs_reverse;
	return x0;
}

// Check if the trajectory_SCP is feasible.
int Planner_MIM::feasible(int order, const double* x0, const Constraint& constraint, const arc* arcs, int length) {
	double* x = new double[order];
	for (int k = 0; k < order; k++) {
		if (!isinf(constraint.M_min[k + 1]) && x0[k] + constraint.EPSILON < constraint.M_min[k + 1]) {
			delete[] x;
			return -(k + 1);
		}
		if (!isinf(constraint.M_max[k + 1]) && x0[k] - constraint.EPSILON > constraint.M_max[k + 1]) {
			delete[] x;
			return (k + 1);
		}
		x[k] = x0[k];
	}
	double* x_old = new double[order];
	double* coeffs = new double[order]();
	double* T = new double[order](); // T[k] = 1/k!
	T[0] = 1;
	double* x_root = new double[order + 1]();
	for (int i = 0; i < length; i++) {
		if (arcs[i].tangent != 0) {
			continue;
		}
		double u = constraint.cal_u(arcs[i]);
		// Check the local minimum or maximum of x[k] along the trajectory_SCP.
		for (int k = 1; k < order; k++) {
			T[k] = T[k - 1] / k;
			if (isinf(constraint.M_max[k]) && isinf(constraint.M_min[k])) {
				continue;
			}
			// Find the root of x[k-1]: $k$th-order polynomial
			for (int j = 0; j < k; j++) {
				coeffs[j] = T[j] * x[k - j - 1];
			}
			coeffs[k] = T[k] * u;
			int num_root = solve_polynomial_real(k, coeffs, x_root);
			for (int l = 0; l < num_root; l++) {
				if (x_root[l] * (x_root[l] - arcs[i].time) < 0) {
					// A local minimum or maximum of x[k] is found.
					double xk = x[k];
					double Tk = 1;
					for (int j = 1; j <= k; j++) {
						Tk *= x_root[l] / j;
						xk += Tk * x[k - j];
					}
					Tk *= x_root[l] / (k + 1);
					xk += Tk * u;
					if (!isinf(constraint.M_min[k + 1]) && xk + constraint.EPSILON < constraint.M_min[k + 1]) {
						delete[] x;
						delete[] x_old;
						delete[] coeffs;
						delete[] x_root;
						return -(k + 1);
					}
					if (!isinf(constraint.M_max[k + 1]) && xk - constraint.EPSILON > constraint.M_max[k + 1]) {
						delete[] x;
						delete[] x_old;
						delete[] coeffs;
						delete[] x_root;
						return k + 1;
					}
				}
			}
		}
		for (int k = 0; k < order; k++) {
			x_old[k] = x[k];
		}
		dynamics_onestep(order, x_old, x, u, arcs[i].time);
		// Check the feasibility of x[k] at the end of the arc.
		for (int k = 0; k < order; k++) {
			if (!isinf(constraint.M_min[k + 1]) && x[k] + constraint.EPSILON < constraint.M_min[k + 1]) {
				delete[] x;
				delete[] x_old;
				delete[] coeffs;
				delete[] x_root;
				return -(k + 1);
			}
			if (!isinf(constraint.M_max[k + 1]) && x[k] - constraint.EPSILON > constraint.M_max[k + 1]) {
				delete[] x;
				delete[] x_old;
				delete[] coeffs;
				delete[] x_root;
				return k + 1;
			}
		}
	}
	delete[] x;
	delete[] x_old;
	delete[] coeffs;
	delete[] x_root;
	return 0;
}

// Determine the direction of 3rd-order trajectory_SCP without position constraint.
int Planner_MIM::direction_3rd_order_no_position(double* x0, double* xf, Constraint& constraint, int direction_2nd) {
	// Step 1. determine direction_2nd
	if (direction_2nd == 0) {
		arc arc_try = plan_1st_order(x0[0], xf[0], constraint);
		double* xs_try = end_point(3, x0, constraint, &arc_try, 1);
		if (norm_inf(xs_try, xf, 2) < constraint.EPSILON) {
			int direction_3rd = xs_try[2] <= xf[2] ? 1 : -1;
			delete[] xs_try;
			return direction_3rd;
		}
		direction_2nd = xs_try[1] < xf[1] ? 1 : -1;
		delete[] xs_try;
	}
	if (direction_2nd < 0) { // x0 is higher than xf. It must be feasible.
		constraint.change_sign(3);
		change_sign_x(x0, 3);
		change_sign_x(xf, 3);
		int direction_3rd = direction_3rd_order_no_position(x0, xf, constraint, 1);
		constraint.change_sign(3);
		change_sign_x(x0, 3);
		change_sign_x(xf, 3);
		return -direction_3rd;
	}
	// x0 is lower than xf
	// Step 2. determine multi_solution
	double multi_solution = false;
	if (xf[0] * x0[0] > constraint.EPSILON) {
		bool change_sign = xf[0] < constraint.EPSILON;
		if (change_sign) {
			constraint.change_sign(2);
			change_sign_x(x0, 2);
			change_sign_x(xf, 2);
		}
		// xf[0]>0
		if (x0[1] - x0[0] * x0[0] / (2.0 * constraint.M_min[0]) < xf[1] - xf[0] * xf[0] / (2.0 * constraint.M_min[0]) && x0[1] - x0[0] * x0[0] / (2.0 * constraint.M_min[0]) > xf[1] - xf[0] * xf[0] / (2.0 * constraint.M_max[0])) {
			if (x0[0] >= xf[0] || x0[1] - x0[0] * x0[0] / (2.0 * constraint.M_max[0]) <= xf[1] - xf[0] * xf[0] / (2.0 * constraint.M_max[0])) {
				multi_solution = true;
			}
		}
		if (change_sign) {
			constraint.change_sign(2);
			change_sign_x(x0, 2);
			change_sign_x(xf, 2);
		}
	}
	// Step 3. Find the proper position
	vector<arc> arcs0 = plan_2nd_order(x0, xf, constraint, false, direction_2nd);
	double* xs0 = end_point(3, x0, constraint, arcs0.data(), arcs0.size());
	double pos0 = xs0[2] - xf[2];
	delete[] xs0;
	if (!multi_solution) {
		return pos0 <= 0 ? 1 : -1;
	}
	// multi_solution = true;
	// Consider 00
	// coeffs[0] + coeffs[1]* T2 + coeffs[2] * T2^2 = 0
	double coeffs_00[3] = { xf[0] * xf[0] - x0[0] * x0[0] + 2.0 * constraint.M_min[0] * (x0[1] - xf[1]),2.0 * xf[0] * (constraint.M_min[0] - constraint.M_max[0]),constraint.M_max[0] * (constraint.M_max[0] - constraint.M_min[0]) };
	static double x_00[2];
	int num_root = solve_polynomial_real(2, coeffs_00, x_00);
	double T21 = min(x_00[0], x_00[1]); // T21 > 0 holds
	double T22 = max(x_00[0], x_00[1]); // T22 > 0 holds
	double T11 = (xf[0] - x0[0] - (constraint.M_max[0] - constraint.M_min[0]) * T21) / constraint.M_min[0]; // T11 > T21 holds
	double T12 = (xf[0] - x0[0] - (constraint.M_max[0] - constraint.M_min[0]) * T22) / constraint.M_min[0]; // T12 > T22 holds
	vector<arc> arcs;
	arcs.push_back(arc(0, false, 0, T11 - T21));
	arcs.push_back(arc(0, true, 0, T21));
	double* xs1 = end_point(3, x0, constraint, arcs.data(), arcs.size());
	double pos1 = xs1[2] - xf[2];
	delete[] xs1;
	arcs.clear();
	if (isinf(constraint.M_min[1]) || T22 <= (xf[0] - constraint.M_min[1]) / constraint.M_max[0]) {
		arcs.push_back(arc(0, false, 0, T12 - T22));
		arcs.push_back(arc(0, true, 0, T22));
	}
	else {
		double t0 = (constraint.M_min[1] - x0[0]) / constraint.M_min[0];
		double t2 = (xf[0] - constraint.M_min[1]) / constraint.M_max[0];
		double t1 = ((xf[1] - xf[0] * t2 + 0.5 * constraint.M_max[0] * t2 * t2) - (x0[1] + x0[0] * t0 + 0.5 * constraint.M_min[0] * t0 * t0)) / constraint.M_min[1];
		arcs.push_back(arc(0, false, 0, t0));
		arcs.push_back(arc(1, false, 0, t1));
		arcs.push_back(arc(0, true, 0, t2));
	}
	double* xs2 = end_point(3, x0, constraint, arcs.data(), arcs.size());
	double pos2 = xs2[2] - xf[2];
	delete[] xs2;
	return ((pos0 <= 0) + (pos1 <= 0) + (pos2 <= 0)) % 2 == 1 ? 1 : -1;
}


// Plan a 3rd-order trajectory_SCP from x0 to xf under constraint. (No position constraint)
// First Step: 00, 010
// Second Step: s00, s010
// If direction is 1/-1, then directly consider the case where x0 is lower/higher than xf.
vector<arc> Planner_MIM::plan_3rd_order_no_position(double* x0, double* xf, Constraint& constraint, int direction) {
	if (direction == 0) {
		direction = direction_3rd_order_no_position(x0, xf, constraint, 0);
	}
	if (direction < 0) { // x0 is higher than xf. It must be feasible.
		constraint.change_sign(3);
		change_sign_x(x0, 3);
		change_sign_x(xf, 3);
		vector<arc> result = plan_3rd_order_no_position(x0, xf, constraint, 1);
		constraint.change_sign(3);
		change_sign_x(x0, 3);
		change_sign_x(xf, 3);
		change_sign_arcs(result.data(), result.size());
		return result;
	}
	// x0 is lower than xf
	// Step 1: Try to enter the manifold.
	vector<arc> result;
	waypoints xs_before;
	if (!isinf(constraint.M_max[2])) {
		double maxV[] = { 0, constraint.M_max[2] };
		result = plan_2nd_order(x0, maxV, constraint, false, 1);
		vector<arc> result2 = plan_2nd_order(maxV, xf, constraint, false, -1);
		waypoints xs1 = end_points(3, x0, constraint, result.data(), result.size());
		waypoints xs2 = start_points(3, xf, constraint, result2.data(), result2.size());
		if (xs1[xs1.length - 1][2] <= xs2[xs2.length - 1][2]) {
			if (xs1[xs1.length - 1][2] < xs2[xs2.length - 1][2]) {
				result.push_back(arc(2, true, 0, (xs2[xs2.length - 1][2] - xs1[xs1.length - 1][2]) / constraint.M_max[2]));
			}
			for (int i = 0; i < result2.size(); i++) {
				result.push_back(result2[i]);
			}
			xs2.clear();
			xs1.clear();
			return result;
		}
		xs1.length -= 2; // +0 cannot exist.
		result.pop_back();
		result[result.size() - 1].time = 0;
		xs_before = xs1;
		xs1.xs = nullptr;
		xs1.length = 0;
		xs1.order = 0;
	}
	else if (!isinf(constraint.M_max[1])) {
		result.push_back(plan_1st_order(x0[0], constraint.M_max[1], constraint));
		result.push_back(arc(1, true, 0, 0));
		xs_before.set(3, 2);
		for (int k = 0; k < 2; k++) {
			xs_before[0][k] = x0[k];
		}
		dynamics_onestep(3, x0, xs_before[1], constraint.M_max[0], result[0].time);
	}
	else {
		result.push_back(arc(0, true, 0, 0));
		xs_before.set(3, 1);
		for (int k = 0; k < 2; k++) {
			xs_before[0][k] = x0[k];
		}
	}
	// Step 2: Find where the trajectory_SCP enters the manifold.
	// printf("%d %d\n", direction_3rd_order_no_position(x0, xf, constraint, 0), direction_3rd_order_no_position(xs_before[1], xf, constraint, 0)); 
	while (result.size() > 1) {
		if (direction_3rd_order_no_position(xs_before[xs_before.length - 1], xf, constraint, 0) > 0) {
			break;
		}
		result.pop_back();
		result[result.size() - 1].time = 0;
		--xs_before.length;
	}
	// xs1[xs1.length-1] enters the manifold along the arc result[result.size()-1].
	// Step 3: Create the whole trajectory_SCP. Consider s00, s010
	// Step 3.1: s00
	arc arc2(0, false, 0);
	arc arc3(0, true, 0);
	//xs_before.length = 1;
	//result.pop_back();
	double* T_s00 = solution_3arc_3rd_order(xs_before[xs_before.length - 1], xf, result[result.size() - 1], arc2, arc3, constraint);
	if (T_s00) {
		result[result.size() - 1].time = T_s00[0] - T_s00[1];
		arc2.time = T_s00[1] - T_s00[2];
		arc3.time = T_s00[2];
		delete[] T_s00;
		result.push_back(arc2);
		result.push_back(arc3);
		if (feasible(2, xs_before[xs_before.length - 1], constraint, result.data() + (result.size() - 3), 3) == 0) {
			return result;
		}
		else {
			result.pop_back();
			result.pop_back();
		}
	}
	else {
		throw std::invalid_argument("No solution for s00.");
		exit(1);
	}
	if (isinf(constraint.M_min[1])) {
		throw std::invalid_argument("No solution for s01.");
		exit(1);
	}
	// Step 3.2: s010
	arc arc4 = plan_1st_order(constraint.M_min[1], xf[0], constraint);
	double* xsf = start_point(3, xf, constraint, &arc4, 1);
	arc3 = arc(1, false, 0);
	double* T_s01 = solution_3arc_3rd_order(xs_before[xs_before.length - 1], xsf, result[result.size() - 1], arc2, arc3, constraint);
	delete[] xsf;
	if (T_s01 != nullptr) {
		result[result.size() - 1].time = T_s01[0] - T_s01[1];
		arc2.time = T_s01[1] - T_s01[2];
		arc3.time = T_s01[2];
		delete[] T_s01;
		result.push_back(arc2);
		result.push_back(arc3);
		if (feasible(2, xs_before[xs_before.length - 1], constraint, result.data() + (result.size() - 3), 3) == 0) {
			result.push_back(arc4);
			return result;
		}
		else {
			result.pop_back();
			result.pop_back();
		}
	}
	throw std::invalid_argument("No solution.");
	exit(1);
}

// Consider tangent marker.
vector<arc> Planner_MIM::tangent_marker_3rd_order(double* x0, Constraint& constraint, bool flag_consider_position, bool forward, bool upper_bound) {
	vector<arc> result;
	// Step 1. Consider x0 -> 00 -> (3,2)
	arc arc1(0, forward != upper_bound, 0);
	arc arc2(0, !arc1.sign, 0);
	double* T_tangent = solution_2arc_tangent_3rd_order(x0, arc1, arc2, constraint, forward);
	if (T_tangent != nullptr) {
		arc1.time = T_tangent[0] - T_tangent[1];
		arc2.time = T_tangent[1];
		delete[] T_tangent;
		result.push_back(arc1);
		result.push_back(arc2);
		if (feasible(2, x0, constraint, result.data(), result.size()) != 0) {
			result.clear();
		}
		else {
			return result;
		}
	}
	// Step 4.1.2: 010(-3,2)*
	arc2 = arc(1, arc1.sign == forward, 0);
	double x1[3];
	x1[0] = arc2.sign ? constraint.M_max[1] : constraint.M_min[1];
	if (!isinf(x1[0])) {
		arc arc3(0, !arc1.sign, 0);
		arc1.time = (x1[0] - x0[0]) / constraint.cal_u(arc1);
		dynamics_onestep(3, x0, x1, constraint.cal_u(arc1), arc1.time);
		T_tangent = solution_2arc_tangent_3rd_order(x1, arc2, arc3, constraint, true);
		if (T_tangent != nullptr) {
			arc2.time = T_tangent[0] - T_tangent[1];
			arc3.time = T_tangent[1];
			delete[] T_tangent;
			result.push_back(arc1);
			result.push_back(arc2);
			result.push_back(arc3);
			if (feasible(2, x0, constraint, result.data(), result.size()) != 0) {
				result.clear();
			}
			else {
				return result;
			}
		}
	}
	result.clear();
	return result;
}

// Plan a 3rd-order trajectory_SCP from x0 to xf under constraint.
// First Step: 00, 010
// Second Step: s00, s010
// If direction is 1/-1, then directly consider the case where x0 is lower/higher than xf.
vector<arc> Planner_MIM::plan_3rd_order(double* x0, double* xf, Constraint& constraint, bool flag_consider_position) {
	vector<arc> result0 = plan_3rd_order_no_position(x0, xf, constraint, 0);
	if ((!flag_consider_position)) {
		return result0;
	}
	int active_constraint = feasible(3, x0, constraint, result0.data(), result0.size());
	if (active_constraint == 0) {
		return result0;
	}
	if (abs(active_constraint) != 3) {
		throw std::invalid_argument("The trajectory_SCP is not feasible in 3rd order.");
		exit(1);
	}
	bool upper_bound = active_constraint == 3;

	// 00(3,2)000
	vector<arc> result10 = tangent_marker_3rd_order(x0, constraint, flag_consider_position, true, upper_bound);
	if (!result10.empty()) {
		double* xs10 = end_point(3, x0, constraint, result10.data(), result10.size());
		vector<arc> result11 = plan_3rd_order_no_position(xs10, xf, constraint, result10.back().sign ? 1 : -1);
		delete[] xs10;
		if (!result11.empty()) {
			result10.push_back(arc(3, upper_bound, 2));
			for (int i = 0; i < result11.size(); i++) {
				result10.push_back(result11[i]);
			}
			if (feasible(3, x0, constraint, result10.data(), result10.size()) == 0) {
				return result10;
			}
		}
	}

	// 000(3,2)00
	vector<arc> result20 = tangent_marker_3rd_order(xf, constraint, flag_consider_position, false, upper_bound);
	if (!result20.empty()) {
		double* xs20 = end_point(3, xf, constraint, result20.data(), result20.size());
		vector<arc> result21 = plan_3rd_order_no_position(x0, xs20, constraint, result20.back().sign ? 1 : -1);
		delete[] xs20;
		if (!result21.empty()) {
			result21.push_back(arc(3, upper_bound, 2));
			for (int i = result20.size() - 1; i >= 0; i--) {
				result20[i].time *= -1;
				result21.push_back(result20[i]);
			}
			if (feasible(3, x0, constraint, result21.data(), result21.size()) == 0) {
				return result21;
			}
		}
	}
	// TODO: more than one tangent markers (small probability)
	return vector <arc>();
}

//// Plan a 3rd-order trajectory_SCP from x0 to xf under constraint.
//// First Step: 00, 010
//// Second Step: s00, s010
//// If direction is 1/-1, then directly consider the case where x0 is lower/higher than xf.
//vector<arc> Planner_MIM::plan_3rd_order(double* x0, double* xf, Constraint& constraint, bool flag_consider_position, int direction) {
//	if (direction == 0) {
//		vector<arc> arc_try = plan_2nd_order(x0, xf, constraint, false, 0);
//		if (arc_try.empty()) {
//			return arc_try;
//		}
//		waypoints xs_try = end_points(3, x0, constraint, arc_try.data(), arc_try.size());
//		if (norm_inf(xs_try[xs_try.length - 1], xf, 3) < constraint.EPSILON) {
//			if (flag_consider_position && (feasible(3, x0, constraint, arc_try.data(), arc_try.size()) != 0)) {
//				arc_try.clear();
//			}
//			return arc_try;
//		}
//		direction = xs_try[xs_try.length - 1][2] < xf[2] ? 1 : -1;
//	}
//	if (direction < 0) { // x0 is higher than xf. It must be feasible.
//		constraint.change_sign(3);
//		change_sign_x(x0, 3);
//		change_sign_x(xf, 3);
//		vector<arc> result = plan_3rd_order(x0, xf, constraint, flag_consider_position, 1);
//		constraint.change_sign(3);
//		change_sign_x(x0, 3);
//		change_sign_x(xf, 3);
//		change_sign_arcs(result.data(), result.size());
//		return result;
//	}
//	// x0 is lower than xf
//	// Step 1: Try to enter the manifold.
//	vector<arc> result;
//	waypoints xs_before;
//	if (!isinf(constraint.M_max[2])) {
//		double maxV[] = { 0, constraint.M_max[2] };
//		result = plan_2nd_order(x0, maxV, constraint, false, 1);
//		vector<arc> result2 = plan_2nd_order(maxV, xf, constraint, false, -1);
//		waypoints xs1 = end_points(3, x0, constraint, result.data(), result.size());
//		waypoints xs2 = start_points(3, xf, constraint, result2.data(), result2.size());
//		if (xs1[xs1.length - 1][2] <= xs2[xs2.length - 1][2]) {
//			if (xs1[xs1.length - 1][2] < xs2[xs2.length - 1][2]) {
//				result.push_back(arc(2, true, 0, (xs2[xs2.length - 1][2] - xs1[xs1.length - 1][2]) / constraint.M_max[2]));
//			}
//			for (int i = 0; i < result2.size(); i++) {
//				result.push_back(result2[i]);
//			}
//			xs2.clear();
//			xs1.clear();
//			return result;
//		}
//		xs1.length -= 2; // +0 cannot exist.
//		result.pop_back();
//		result[result.size() - 1].time = 0;
//		xs_before = xs1;
//		xs1.xs = nullptr;
//		xs1.length = 0;
//		xs1.order = 0;
//	}
//	else if (!isinf(constraint.M_max[1])) {
//		result.push_back(plan_1st_order(x0[0], constraint.M_max[1], constraint));
//		result.push_back(arc(1, true, 0, 0));
//		xs_before.set(3, 2);
//		for (int k = 0; k < 2; k++) {
//			xs_before[0][k] = x0[k];
//		}
//		dynamics_onestep(3, x0, xs_before[1], constraint.M_max[0], result[0].time);
//	}
//	else {
//		result.push_back(arc(0, true, 0, 0));
//		xs_before.set(3, 1);
//		for (int k = 0; k < 2; k++) {
//			xs_before[0][k] = x0[k];
//		}
//	}
//	// Step 2: Find where the trajectory_SCP enters the manifold.
//	bool flag_succeed_noposition = false;
//	while (result.size() > 1) {
//		double propos = proper_position(3, xs_before[xs_before.length - 1], xf, constraint);
//		if (xs_before[xs_before.length - 1][2] <= xf[2] - propos) {
//			break;
//		}
//		result.pop_back();
//		result[result.size() - 1].time = 0;
//		--xs_before.length;
//	}
//	// xs1[xs1.length-1] enters the manifold along the arc result[result.size()-1].
//	// Step 3: Create the whole trajectory_SCP. Consider s00, s010
//	// Step 3.1: s00
//	arc arc2(0, false, 0);
//	arc arc3(0, true, 0);
//	//xs_before.length = 1;
//	//result.pop_back();
//	double* T_s00 = solution_3arc_3rd_order(xs_before[xs_before.length - 1], xf, result[result.size() - 1], arc2, arc3, constraint);
//	if (T_s00) {
//		result[result.size() - 1].time = T_s00[0] - T_s00[1];
//		arc2.time = T_s00[1] - T_s00[2];
//		arc3.time = T_s00[2];
//		delete[] T_s00;
//		result.push_back(arc2);
//		result.push_back(arc3);
//		if (feasible(2, xs_before[xs_before.length - 1], constraint, result.data() + (result.size() - 3), 3) == 0) {
//			flag_succeed_noposition = true;
//		}
//		else {
//			result.pop_back();
//			result.pop_back();
//		}
//	}
//	else {
//		throw std::invalid_argument("No solution for s00.");
//		exit(1);
//	}
//	if (!flag_succeed_noposition) {
//		if (isinf(constraint.M_min[1])) {
//			throw std::invalid_argument("No solution for s01.");
//			exit(1);
//		}
//		// Step 3.2: s010
//		arc arc4 = plan_1st_order(constraint.M_min[1], xf[0], constraint);
//		waypoints xsf = start_points(3, xf, constraint, &arc4, 1);
//		arc3 = arc(1, false, 0);
//		double* T_s01 = solution_3arc_3rd_order(xs_before[xs_before.length - 1], xsf[xsf.length - 1], result[result.size() - 1], arc2, arc3, constraint);
//		if (T_s01 != nullptr) {
//			result[result.size() - 1].time = T_s01[0] - T_s01[1];
//			arc2.time = T_s01[1] - T_s01[2];
//			arc3.time = T_s01[2];
//			delete[] T_s01;
//			result.push_back(arc2);
//			result.push_back(arc3);
//			if (feasible(2, xs_before[xs_before.length - 1], constraint, result.data() + (result.size() - 3), 3) == 0) {
//				flag_succeed_noposition = true;
//				result.push_back(arc4);
//			}
//			else {
//				result.pop_back();
//				result.pop_back();
//			}
//		}
//	}
//	if (!flag_succeed_noposition) {
//		throw std::invalid_argument("No solution.");
//		exit(1);
//	}
//	if (feasible(3, x0, constraint, result.data(), result.size()) == 0) {
//		return result;
//	}
//	// Step 4: Consider tangent markers.
//	// Only -3 can support tangent markers, i.e., (-3,2).
//	result.clear();
//	// Step 4.1: like 00(-3,2)000
//	// Step 4.1.1: 00(-3,2)*
//	arc arc1(0, true, 0);
//	arc2 = arc(0, false, 0);
//	double* T_tangent = solution_2arc_tangent_3rd_order(x0, arc1, arc2, constraint, true);
//	if (T_tangent != nullptr) {
//		arc1.time = T_tangent[0] - T_tangent[1];
//		arc2.time = T_tangent[1];
//		delete[] T_tangent;
//		result.push_back(arc1);
//		result.push_back(arc2);
//		if (feasible(2, x0, constraint, result.data(), result.size())) {
//			waypoints xs1 = end_points(3, x0, constraint, result.data(), result.size());
//			result.push_back(arc(3, true, 2));
//			vector<arc> result2 = plan_3rd_order(xs1[2], xf, constraint, false, -1);
//			for (int i = 0; i < result2.size(); ++i) {
//				result.push_back(result2[i]);
//			}
//			return result;
//		}
//		else {
//			result.clear();
//		}
//	}
//	// Step 4.1.2: 010(-3,2)*
//	arc1 = plan_1st_order(x0[0], constraint.M_max[1], constraint);
//	arc2 = arc(1, true, 0);
//	arc3 = arc(0, false, 0);
//	waypoints x0now_pre = end_points(3, x0, constraint, &arc1, 1);
//	T_tangent = solution_2arc_tangent_3rd_order(x0now_pre[1], arc2, arc3, constraint, true);
//	if (T_tangent != nullptr) {
//		arc2.time = T_tangent[0] - T_tangent[1];
//		arc3.time = T_tangent[1];
//		delete[] T_tangent;
//		result.push_back(arc1);
//		result.push_back(arc2);
//		result.push_back(arc3);
//		if (feasible(2, x0, constraint, result.data(), result.size())) {
//			waypoints xs1 = end_points(3, x0, constraint, result.data(), result.size());
//			result.push_back(arc(3, true, 2));
//			vector<arc> result2 = plan_3rd_order(xs1[2], xf, constraint, false, -1);
//			for (int i = 0; i < result2.size(); ++i) {
//				result.push_back(result2[i]);
//			}
//		}
//		return result;
//	}
//	// Step 4.2: like 000(-3,2)00
//	// Step 4.2.1: *(-3,2)00
//	arc1 = arc(0, false, 0);
//	arc2 = arc(0, true, 0);
//	T_tangent = solution_2arc_tangent_3rd_order(xf, arc1, arc2, constraint, false);
//	if (T_tangent != nullptr) {
//		arc1.time = T_tangent[1] - T_tangent[0];
//		arc2.time = -T_tangent[1];
//		delete[] T_tangent;
//		arc arcs2[] = { arc2, arc1 };
//		waypoints xs1 = start_points(3, xf, constraint, arcs2, 2);
//		if (feasible(2, xs1[xs1.length - 1], constraint, result.data(), result.size())) {
//			vector<arc> result2 = plan_3rd_order(x0, xs1[2], constraint, false, -1);
//			for (int i = 0; i < result2.size(); ++i) {
//				result.push_back(result2[i]);
//			}
//			result.push_back(arc(3, true, 2));
//			result.push_back(arc2);
//			result.push_back(arc1);
//			return result;
//		}
//	}
//	// Step 4.2.1: *(-3,2)010
//	arc1 = plan_1st_order(constraint.M_max[1], xf[0], constraint);
//	arc2 = arc(1, true, 0);
//	arc3 = arc(0, true, 0);
//	waypoints x0now_post = start_points(3, xf, constraint, &arc1, 1);
//	T_tangent = solution_2arc_tangent_3rd_order(x0now_pre[1], arc2, arc3, constraint, false);
//	if (T_tangent != nullptr) {
//		arc2.time = T_tangent[1] - T_tangent[0];
//		arc3.time = -T_tangent[1];
//		delete[] T_tangent;
//		arc arcs2[] = { arc3, arc2, arc1 };
//		waypoints xs1 = start_points(3, xf, constraint, arcs2, 3);
//		if (feasible(2, xs1[xs1.length - 1], constraint, result.data(), result.size())) {
//			vector<arc> result2 = plan_3rd_order(x0, xs1[3], constraint, false, -1);
//			for (int i = 0; i < result2.size(); ++i) {
//				result.push_back(result2[i]);
//			}
//			result.push_back(arc(3, true, 2));
//			result.push_back(arc3);
//			result.push_back(arc2);
//			result.push_back(arc1);
//			return result;
//		}
//	}
//	result.clear();
//	return result;
//}

// Calculate the infinity norm of x-y.
double Planner_MIM::norm_inf(double* x, double* y, int n) {
	double dmax = 0;
	for (int i = 0; i < n; i++) {
		dmax = ::std::max(dmax, ::std::abs(x[i] - y[i]));
	}
	return dmax;
}

// Calculate the proper position of the trajectory_SCP.
double Planner_MIM::proper_position(int order, double* x0, double* xf, Constraint& constraint) {
	vector<arc> arcs = plan(order - 1, x0, xf, constraint, false);
	double* x0_ = new double[order];
	for (int k = 0; k < order - 1; k++) {
		x0_[k] = x0[k];
	}
	x0_[order - 1] = 0;
	double* xs = end_point(order, x0_, constraint, arcs.data(), arcs.size());
	double pos = xs[order - 1];
	delete[] x0_;
	delete[] xs;
	return pos;
}

// Solve 3arcs in 3rd-order problems based on Groebner basis.
double* Planner_MIM::solution_3arc_3rd_order(double* x0, double* xf, const arc& arc1, const arc& arc2, const arc& arc3, const Constraint& constraint) {
	double u1 = constraint.cal_u(arc1);
	double u2 = constraint.cal_u(arc2);
	double u3 = constraint.cal_u(arc3);
	u3 -= u2;
	u2 -= u1;
	double x1 = x0[0];
	double x2 = x0[1];
	double y1 = xf[0] - x0[0];
	double y2 = xf[1] - x0[1];
	double y3 = xf[2] - x0[2];

	if (::std::abs(u2 - u3) > 1e-6) {
		double c6 = u1 * u1 * (u1 + u2) * (u1 + u3) * (u1 + u2 + u3) * (u1 + u2 + u3);
		double c5 = 6 * u1 * (u1 + u2) * (u1 + u3) * (u2 * x1 + u3 * x1 - u1 * y1) * (u1 + u2 + u3);
		double c4 = (15 * y1 * y1 - 6 * u2 * y2 - 6 * u3 * y2) * u1 * u1 * u1 * u1 + (18 * u2 * y1 * y1 + 18 * u3 * y1 * y1 - 12 * u2 * u2 * y2 - 12 * u3 * u3 * y2 + 24 * u2 * u3 * x2 - 12 * u2 * u3 * y2 - 24 * u2 * x1 * y1 - 24 * u3 * x1 * y1) * u1 * u1 * u1 + (12 * u2 * u2 * x1 * x1 + 12 * u3 * u3 * x1 * x1 + 3 * u2 * u2 * y1 * y1 + 3 * u3 * u3 * y1 * y1 - 6 * u2 * u2 * u2 * y2 - 6 * u3 * u3 * u3 * y2 + 12 * u2 * u3 * x1 * x1 + 36 * u2 * u3 * u3 * x2 + 36 * u2 * u2 * u3 * x2 + 15 * u2 * u3 * y1 * y1 - 6 * u2 * u3 * u3 * y2 - 6 * u2 * u2 * u3 * y2 - 24 * u2 * u2 * x1 * y1 - 24 * u3 * u3 * x1 * y1 - 60 * u2 * u3 * x1 * y1) * u1 * u1 + (12 * x2 * u2 * u2 * u2 * u3 + 12 * u2 * u2 * u2 * x1 * x1 + 24 * x2 * u2 * u2 * u3 * u3 + 24 * u2 * u2 * u3 * x1 * x1 - 30 * y1 * u2 * u2 * u3 * x1 + 12 * x2 * u2 * u3 * u3 * u3 + 24 * u2 * u3 * u3 * x1 * x1 - 30 * y1 * u2 * u3 * u3 * x1 + 12 * u3 * u3 * u3 * x1 * x1) * u1 + 9 * u2 * u2 * u2 * u3 * x1 * x1 + 18 * u2 * u2 * u3 * u3 * x1 * x1 + 9 * u2 * u3 * u3 * u3 * x1 * x1;
		double c3 = (24 * u2 * y1 * y2 - 24 * u2 * u3 * y3 - 20 * y1 * y1 * y1 + 24 * u3 * y1 * y2) * u1 * u1 * u1 + (36 * u2 * x1 * y1 * y1 - 12 * u3 * y1 * y1 * y1 - 36 * u2 * u3 * u3 * y3 - 36 * u2 * u2 * u3 * y3 - 12 * u2 * y1 * y1 * y1 + 36 * u3 * x1 * y1 * y1 - 24 * u2 * u2 * x1 * y2 - 24 * u3 * u3 * x1 * y2 + 24 * u2 * u2 * y1 * y2 + 24 * u3 * u3 * y1 * y2 - 24 * u2 * u3 * x1 * y2 - 72 * u2 * u3 * x2 * y1 + 24 * u2 * u3 * y1 * y2) * u1 * u1 + (12 * u2 * u2 * x1 * y1 * y1 - 12 * u2 * u3 * u3 * u3 * y3 - 12 * u2 * u2 * u2 * u3 * y3 - 24 * u2 * u2 * u2 * x1 * y2 - 24 * u3 * u3 * u3 * x1 * y2 - 24 * u2 * u2 * u3 * u3 * y3 - 4 * u2 * u3 * y1 * y1 * y1 - 24 * u2 * u2 * x1 * x1 * y1 + 12 * u3 * u3 * x1 * y1 * y1 - 24 * u3 * u3 * x1 * x1 * y1 + 72 * u2 * u3 * u3 * x1 * x2 + 72 * u2 * u2 * u3 * x1 * x2 + 48 * u2 * u3 * x1 * y1 * y1 - 24 * u2 * u3 * x1 * x1 * y1 - 12 * u2 * u3 * u3 * x1 * y2 - 36 * u2 * u3 * u3 * x2 * y1 - 12 * u2 * u2 * u3 * x1 * y2 - 36 * u2 * u2 * u3 * x2 * y1 + 12 * u2 * u3 * u3 * y1 * y2 + 12 * u2 * u2 * u3 * y1 * y2) * u1 + 36 * x2 * u2 * u2 * u2 * u3 * x1 + 8 * u2 * u2 * u2 * x1 * x1 * x1 + 72 * x2 * u2 * u2 * u3 * u3 * x1 - 8 * u2 * u2 * u3 * x1 * x1 * x1 - 36 * y1 * u2 * u2 * u3 * x1 * x1 + 36 * x2 * u2 * u3 * u3 * u3 * x1 - 8 * u2 * u3 * u3 * x1 * x1 * x1 - 36 * y1 * u2 * u3 * u3 * x1 * x1 + 8 * u3 * u3 * u3 * x1 * x1 * x1;
		double c2 = (12 * u2 * u2 * y2 * y2 + 72 * y3 * u2 * u3 * y1 + 12 * u2 * u3 * y2 * y2 - 36 * u2 * y1 * y1 * y2 + 12 * u3 * u3 * y2 * y2 - 36 * u3 * y1 * y1 * y2 + 15 * y1 * y1 * y1 * y1) * u1 * u1 + (12 * u2 * u2 * u2 * y2 * y2 + 36 * y3 * u2 * u2 * u3 * y1 - 12 * u2 * u2 * u3 * y2 * y2 - 72 * x2 * u2 * u2 * u3 * y2 - 72 * x1 * y3 * u2 * u2 * u3 - 12 * u2 * u2 * y1 * y1 * y2 + 48 * x1 * u2 * u2 * y1 * y2 + 36 * y3 * u2 * u3 * u3 * y1 - 12 * u2 * u3 * u3 * y2 * y2 - 72 * x2 * u2 * u3 * u3 * y2 - 72 * x1 * y3 * u2 * u3 * u3 - 12 * u2 * u3 * y1 * y1 * y2 + 72 * x2 * u2 * u3 * y1 * y1 + 48 * x1 * u2 * u3 * y1 * y2 + 3 * u2 * y1 * y1 * y1 * y1 - 24 * x1 * u2 * y1 * y1 * y1 + 12 * u3 * u3 * u3 * y2 * y2 - 12 * u3 * u3 * y1 * y1 * y2 + 48 * x1 * u3 * u3 * y1 * y2 + 3 * u3 * y1 * y1 * y1 * y1 - 24 * x1 * u3 * y1 * y1 * y1) * u1 - 36 * y3 * u2 * u2 * u2 * u3 * x1 + 36 * u2 * u2 * u2 * u3 * x2 * x2 - 24 * y2 * u2 * u2 * u2 * x1 * x1 - 72 * y3 * u2 * u2 * u3 * u3 * x1 + 72 * u2 * u2 * u3 * u3 * x2 * x2 + 24 * y2 * u2 * u2 * u3 * x1 * x1 - 72 * u2 * u2 * u3 * x1 * x2 * y1 + 36 * y2 * u2 * u2 * u3 * x1 * y1 + 12 * u2 * u2 * x1 * x1 * y1 * y1 - 36 * y3 * u2 * u3 * u3 * u3 * x1 + 36 * u2 * u3 * u3 * u3 * x2 * x2 + 24 * y2 * u2 * u3 * u3 * x1 * x1 - 72 * u2 * u3 * u3 * x1 * x2 * y1 + 36 * y2 * u2 * u3 * u3 * x1 * y1 + 12 * u2 * u3 * x1 * x1 * y1 * y1 - 12 * u2 * u3 * x1 * y1 * y1 * y1 - 24 * y2 * u3 * u3 * u3 * x1 * x1 + 12 * u3 * u3 * x1 * x1 * y1 * y1;
		double c1 = (72 * y3 * u2 * u2 * u3 * y2 - 24 * u2 * u2 * y1 * y2 * y2 + 72 * y3 * u2 * u3 * u3 * y2 - 72 * y3 * u2 * u3 * y1 * y1 - 24 * u2 * u3 * y1 * y2 * y2 + 24 * u2 * y1 * y1 * y1 * y2 - 24 * u3 * u3 * y1 * y2 * y2 + 24 * u3 * y1 * y1 * y1 * y2 - 6 * y1 * y1 * y1 * y1 * y1) * u1 - 72 * x2 * y3 * u2 * u2 * u2 * u3 + 24 * x1 * u2 * u2 * u2 * y2 * y2 - 144 * x2 * y3 * u2 * u2 * u3 * u3 + 72 * x2 * u2 * u2 * u3 * y1 * y2 + 72 * x1 * y3 * u2 * u2 * u3 * y1 - 24 * x1 * u2 * u2 * u3 * y2 * y2 - 24 * x1 * u2 * u2 * y1 * y1 * y2 - 72 * x2 * y3 * u2 * u3 * u3 * u3 + 72 * x2 * u2 * u3 * u3 * y1 * y2 + 72 * x1 * y3 * u2 * u3 * u3 * y1 - 24 * x1 * u2 * u3 * u3 * y2 * y2 - 24 * x2 * u2 * u3 * y1 * y1 * y1 - 24 * x1 * u2 * u3 * y1 * y1 * y2 + 6 * x1 * u2 * y1 * y1 * y1 * y1 + 24 * x1 * u3 * u3 * u3 * y2 * y2 - 24 * x1 * u3 * u3 * y1 * y1 * y2 + 6 * x1 * u3 * y1 * y1 * y1 * y1;
		double c0 = 36 * u2 * u2 * u2 * u3 * y3 * y3 - 8 * u2 * u2 * u2 * y2 * y2 * y2 + 72 * u2 * u2 * u3 * u3 * y3 * y3 - 72 * u2 * u2 * u3 * y1 * y2 * y3 + 8 * u2 * u2 * u3 * y2 * y2 * y2 + 12 * u2 * u2 * y1 * y1 * y2 * y2 + 36 * u2 * u3 * u3 * u3 * y3 * y3 - 72 * u2 * u3 * u3 * y1 * y2 * y3 + 8 * u2 * u3 * u3 * y2 * y2 * y2 + 24 * u2 * u3 * y1 * y1 * y1 * y3 + 12 * u2 * u3 * y1 * y1 * y2 * y2 - 6 * u2 * y1 * y1 * y1 * y1 * y2 - 8 * u3 * u3 * u3 * y2 * y2 * y2 + 12 * u3 * u3 * y1 * y1 * y2 * y2 - 6 * u3 * y1 * y1 * y1 * y1 * y2 + y1 * y1 * y1 * y1 * y1 * y1;

		double coeffs[7] = { c0, c1, c2, c3, c4, c5, c6 };
		static double x[6];
		int num_root = solve_polynomial_real(6, coeffs, x);
		//printf("ours[%d-order]: ", num_root);
		//for (int i = 0; i < num_root; ++i) {
		//	printf("%f ", x[i]);
		//}
		//printf("\n");

		//static double z[12];
		//int n = 7;
		//while (fabs(coeffs[--n]) < 1e-6);
		//gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(n + 1);
		//gsl_poly_complex_solve(coeffs, n + 1, w, z);
		//gsl_poly_complex_workspace_free(w);
		//printf("gsl[%d-order]: ", n);
		//for (int i = 0; i < n; ++i) {
		//	if (fabs(z[2 * i + 1]) < 1e-6) {
		//		printf("%f ", z[2 * i]);
		//	}
		//}
		//printf("\n");

		bool flag_succeed = false;
		double* T = new double[3]();
		for (int i = 0; i < num_root; ++i) {
			if (x[i] > 0) {
				double T1 = x[i];
				double T2_under = -u2 * (-T1 * T1 * T1 * u1 * u1 * u1 * u2 + 3 * T1 * T1 * T1 * u1 * u1 * u1 * u3 - T1 * T1 * T1 * u1 * u1 * u2 * u2 + 3 * T1 * T1 * T1 * u1 * u1 * u2 * u3 + 4 * T1 * T1 * T1 * u1 * u1 * u3 * u3 + T1 * T1 * T1 * u1 * u2 * u2 * u3 + 2 * T1 * T1 * T1 * u1 * u2 * u3 * u3 + T1 * T1 * T1 * u1 * u3 * u3 * u3 + 3 * T1 * T1 * u1 * u1 * u2 * y1 - 9 * T1 * T1 * u1 * u1 * u3 * y1 + T1 * T1 * u1 * u2 * u2 * y1 - 2 * x1 * T1 * T1 * u1 * u2 * u2 - 3 * T1 * T1 * u1 * u2 * u3 * y1 + 6 * x1 * T1 * T1 * u1 * u2 * u3 - 4 * T1 * T1 * u1 * u3 * u3 * y1 + 8 * x1 * T1 * T1 * u1 * u3 * u3 + 3 * x1 * T1 * T1 * u2 * u2 * u3 + 6 * x1 * T1 * T1 * u2 * u3 * u3 + 3 * x1 * T1 * T1 * u3 * u3 * u3 + 2 * y2 * T1 * u1 * u2 * u2 - 6 * y2 * T1 * u1 * u2 * u3 - 3 * T1 * u1 * u2 * y1 * y1 - 8 * y2 * T1 * u1 * u3 * u3 + 9 * T1 * u1 * u3 * y1 * y1 + 6 * x2 * T1 * u2 * u2 * u3 + 2 * x1 * T1 * u2 * u2 * y1 + 12 * x2 * T1 * u2 * u3 * u3 - 6 * x1 * T1 * u2 * u3 * y1 + 6 * x2 * T1 * u3 * u3 * u3 - 8 * x1 * T1 * u3 * u3 * y1 - 6 * y3 * u2 * u2 * u3 - 2 * y2 * u2 * u2 * y1 - 12 * y3 * u2 * u3 * u3 + 6 * y2 * u2 * u3 * y1 + u2 * y1 * y1 * y1 - 6 * y3 * u3 * u3 * u3 + 8 * y2 * u3 * u3 * y1 - 3 * u3 * y1 * y1 * y1);
				double T2_over = (u2 - u3) * (T1 * T1 * u1 * u1 + u3 * T1 * T1 * u1 - 2 * T1 * u1 * y1 + 2 * u3 * x1 * T1 + y1 * y1 - 2 * u3 * y2) * (y1 * y1 - 2 * u3 * y2 - 2 * u2 * y2 + T1 * T1 * u1 * u1 + 2 * T1 * u2 * x1 + 2 * T1 * u3 * x1 - 2 * T1 * u1 * y1 + T1 * T1 * u1 * u2 + T1 * T1 * u1 * u3);
				// T2_under * T2 + T2_over == 0
				double T2 = -T2_over / T2_under;
				double T3 = (u1 * T1 + u2 * T2 - y1) / (-u3);
				if (0 <= T3 && T3 <= T2 && T2 <= T1) {
					if (flag_succeed && T1 <= T[0]) {
						continue;
						// throw std::invalid_argument("Multiple solutions.");
						// exit(1); // TODO
					}
					flag_succeed = true;
					T[0] = T1;
					T[1] = T2;
					T[2] = T3;
					//break;
				}
			}
		}
		if (flag_succeed) {
			return T;
		}
	}
	else {
		double c3 = u1 * (u1 * u1 + 3 * u1 * u2 + 2 * u2 * u2);
		double c2 = (-3 * y1) * u1 * u1 + (6 * u2 * x1 - 3 * u2 * y1) * u1 + 6 * u2 * u2 * x1;
		double c1 = (3 * y1 * y1 - 6 * u2 * y2) * u1 + 12 * x2 * u2 * u2 - 6 * x1 * y1 * u2;
		double c0 = -12 * y3 * u2 * u2 + 6 * y2 * u2 * y1 - y1 * y1 * y1;

		double coeffs[4] = { c0, c1, c2, c3 };
		static double x[3];
		int num_root = solve_polynomial_real(3, coeffs, x);
		//int n = 4;
		//while (fabs(coeffs[--n]) < 1e-6);
		//gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(n + 1);
		//gsl_poly_complex_solve(coeffs, n + 1, w, z);
		//gsl_poly_complex_workspace_free(w);
		bool flag_succeed = false;
		double* T = new double[3]();
		for (int i = 0; i < num_root; ++i) {
			if (x[i] > 0) {
				double T1 = x[i];
				double D2 = 2 * u2 * u2;
				double D1 = 2 * T1 * u1 * u2 - 2 * u2 * y1;
				double D0 = T1 * T1 * u1 * u1 + u2 * T1 * T1 * u1 - 2 * T1 * u1 * y1 + 2 * u2 * x1 * T1 + y1 * y1 - 2 * u2 * y2;
				double Delta = D1 * D1 - 4 * D2 * D0;
				if (Delta < 0) {
					continue;
				}
				double a = ::std::sqrt(float(1.0));
				double T2 = (-D1 + ::std::sqrt(Delta)) / (2 * D2);
				double T3 = (u1 * T1 + u2 * T2 - y1) / (-u2);
				if (0 <= T3 && T3 <= T2 && T2 <= T1) {
					if (flag_succeed) {
						throw "Multiple solutions.";
						exit(1); // TODO
					}
					flag_succeed = true;
					T[0] = T1;
					T[1] = T2;
					T[2] = T3;
					//break;
				}
				T2 = (-D1 - sqrt(Delta)) / (2 * D2);
				T3 = (u1 * T1 + u2 * T2 - y1) / (-u2);
				if (0 <= T3 && T3 <= T2 && T2 <= T1) {
					if (flag_succeed) {
						throw "Multiple solutions.";
						exit(1); // TODO
					}
					flag_succeed = true;
					T[0] = T1;
					T[1] = T2;
					T[2] = T3;
					//break;
				}
			}
		}
		if (flag_succeed) {
			return T;
		}
	}
	//throw "No solutions.";
	//exit(1); // TODO
	return nullptr;
}



//// Solve 3arcs in 3rd-order problems based on Groebner basis.
//double* Planner_MIM::solution_3arc_3rd_order(double* x0, double* xf, const arc& arc1, const arc& arc2, const arc& arc3, const Constraint& constraint) {
//	double u1 = constraint.cal_u(arc1);
//	double u2 = constraint.cal_u(arc2);
//	double u3 = constraint.cal_u(arc3);
//	u3 -= u2;
//	u2 -= u1;
//	double x1 = x0[0];
//	double x2 = x0[1];
//	double y1 = xf[0] - x0[0];
//	double y2 = xf[1] - x0[1];
//	double y3 = xf[2] - x0[2];
//
//	// plus
//	double u12 = u1 + u2;
//	double u13 = u1 + u3;
//	double u23 = u2 + u3;
//	double u123 = u12 + u3;
//	double u1m1 = u1 * u1;
//	double u1m3 = u1 * u3;
//	double u2m2 = u2 * u2;
//	double u2m3 = u2 * u3;
//	double u2m2m2 = u2m2 * u2;
//	double u3m3 = u3 * u3;
//	double u3m3m3 = u3m3 * u3;
//	double x1m1 = x1 * x1;
//	double x1m1m1 = x1m1 * x1;
//	double x2m2 = x2 * x2;
//	double y1m1 = y1 * y1;
//	double y1m1m1 = y1m1 * y1;
//	double y1m2 = y1 * y2;
//	double y2m2 = y2 * y2;
//	double y2m2m2 = y2m2 * y2;
//	double y1m3 = y1 * y3;
//	double u3mx1 = u3 * x1;
//	double u3mx2 = u3 * x2;
//	double u3my2 = u3 * y2;
//	double x1my1 = x1 * y1;
//
//
//	if (fabs(u2 - u3) > 1e-6) {
//		double c6 = u1m1 * u12 * u13 * u123 * u123;
//		double c5 = 6 * u1 * u12 * u13 * (u23 * x1 - u1 * y1) * u123;
//		double c4 = -3 * (-3 * u2m3 * u23 * u23 * x1m1 - 2 * u1 * u23 * (2 * u3m3 * x1m1 + 2 * u2m2 * (x1m1 + u3mx2) + u2m3 * (2 * x1m1 + 2 * u3mx2 - 5 * x1my1)) + u1m1 * u1m1 * (-5 * y1m1 + 2 * u23 * y2) + u1m1 * u1 * (4 * u2m2 * y2 + 2 * u3 * (4 * x1my1 - 3 * y1m1 + 2 * u3my2) + u2 * (-8 * u3mx2 + 8 * x1my1 - 6 * y1m1 + 4 * u3my2)) + u1m1 * (2 * u2m2m2 * y2 - u2m2 * (4 * x1m1 + 12 * u3mx2 - 8 * x1my1 + y1m1 - 2 * u3my2) + u3m3 * (-4 * x1m1 + 8 * x1my1 - y1m1 + 2 * u3my2) + u2m3 * (-4 * x1m1 + 20 * x1my1 - 5 * y1m1 + 2 * u3 * (-6 * x2 + y2))));
//		double c3 = -4 * (-((u3mx1 - u1 * y1) * (5 * u1m1 * y1m1 + u1m3 * (-4 * x1my1 + 3 * y1m1 - 6 * u1 * y2) + 2 * u3m3 * (x1m1 - 3 * u1 * y2))) + u2m2m2 * (-2 * x1m1m1 - 9 * u3mx1 * x2 + 6 * u1 * x1 * y2 + 3 * u1m3 * y3) + u2m2 * (3 * u1 * (2 * x1m1 * y1 - x1 * y1m1 + 2 * u1 * x1 * y2 - 2 * u1 * y1m2) + 6 * u3m3 * (-3 * x1 * x2 + u1 * y3) + u3 * (2 * x1m1m1 + 9 * x1m1 * y1 + 3 * u1 * x1 * (-6 * x2 + y2) + 3 * u1 * (3 * x2 * y1 - y1m2 + 3 * u1 * y3))) + u2 * (3 * u1m1 * y1 * (-3 * x1my1 + y1m1 - 2 * u1 * y2) + u3m3m3 * (-9 * x1 * x2 + 3 * u1 * y3) + u1m3 * (6 * x1m1 * y1 + 18 * u1 * x2 * y1 + y1m1m1 - 6 * u1 * y1m2 + 6 * x1 * (-2 * y1m1 + u1 * y2) + 6 * u1m1 * y3) + u3m3 * (2 * x1m1m1 + 9 * x1m1 * y1 + 3 * u1 * x1 * (-6 * x2 + y2) + 3 * u1 * (3 * x2 * y1 - y1m2 + 3 * u1 * y3))));
//		double c2 = 3 * ((y1m1 - 2 * u3my2) * (5 * u1m1 * y1m1 + u3m3 * (4 * x1m1 - 2 * u1 * y2) + u1m3 * (-8 * x1my1 + y1m1 - 2 * u1 * y2)) + 4 * u2m2m2 * (y2 * (-2 * x1m1 + u1 * y2) + 3 * u3 * (x2m2 - x1 * y3)) + 4 * u2m2 * (x1m1 * y1m1 + 4 * u1 * x1my1 * y2 + u1 * y2 * (-y1m1 + u1 * y2) + 6 * u3m3 * (x2m2 - x1 * y3) + u3 * (2 * x1m1 * y2 + x1 * (-6 * x2 * y1 + 3 * y1m2 - 6 * u1 * y3) - u1 * (6 * x2 * y2 + y2m2 - 3 * y1m3))) + u2 * (u1 * y1m1 * (-8 * x1my1 + y1m1 - 12 * u1 * y2) + 12 * u3m3m3 * (x2m2 - x1 * y3) + 4 * u3m3 * (2 * x1m1 * y2 + x1 * (-6 * x2 * y1 + 3 * y1m2 - 6 * u1 * y3) - u1 * (6 * x2 * y2 + y2m2 - 3 * y1m3)) + 4 * u3 * (x1m1 * y1m1 - x1 * (y1m1m1 - 4 * u1 * y1m2) + u1 * (6 * x2 * y1m1 - y1m1 * y2 + u1 * y2m2 + 6 * u1 * y1m3))));
//		double c1 = 6 * ((u3mx1 - u1 * y1) * (y1m1 - 2 * u3my2) * (y1m1 - 2 * u3my2) + 4 * u2m2m2 * (x1 * y2m2 - 3 * u3mx2 * y3) - 4 * u2m2 * (y1m2 * (x1my1 + u1 * y2) + 6 * u3m3 * x2 * y3 + u3 * (-3 * x2 * y1m2 + x1 * y2m2 - 3 * x1 * y1m3 - 3 * u1 * y2 * y3)) + u2 * (y1m1m1 * (x1my1 + 4 * u1 * y2) - 12 * u3m3m3 * x2 * y3 - 4 * u3 * y1 * (x2 * y1m1 + x1my1 * y2 + u1 * y2m2 + 3 * u1 * y1m3) + 4 * u3m3 * (3 * x2 * y1m2 - x1 * y2m2 + 3 * x1 * y1m3 + 3 * u1 * y2 * y3)));
//		double c0 = y1m1m1 * y1m1m1 - 6 * u23 * y1m1 * y1m1 * y2 + 12 * (u2m2 + u2m3 + u3m3) * y1m1 * y2m2 + 24 * u2m3 * y1m1m1 * y3 - 72 * u2m3 * u23 * y1m2 * y3 + 4 * u23 * (-2 * u3m3 * y2m2m2 + u2m2 * (-2 * y2m2m2 + 9 * u3 * y3 * y3) + u2m3 * (4 * y2m2m2 + 9 * u3 * y3 * y3));
//
//		double coeffs[7] = { c0, c1, c2, c3, c4, c5, c6 };
//		static double x[6];
//		int num_root = solve_polynomial_real(6, coeffs, x);
//		bool flag_succeed = false;
//		double* T = new double[3]();
//		for (int i = 0; i < num_root; ++i) {
//			if (x[2 * i] > 0) {
//				double T1 = x[2 * i];
//				double T1mT1 = T1 * T1;
//				double T2_under = -u2 * (T1mT1 * T1 * u1 * (-(u1m1 * (u2 - 3 * u3)) + u3 * u23 * u23 + u1 * (-u2m2 + 3 * u2m3 + 4 * u3m3)) + T1mT1 * (3 * u3 * u23 * u23 * x1 - u1 * (u2m2 - 3 * u2m3 - 4 * u3m3) * (2 * x1 - y1) + 3 * u1m1 * (u2 - 3 * u3) * y1) + T1 * (2 * u2m2 * (3 * u3mx2 + x1my1 + u1 * y2) + u3 * (6 * u3m3 * x2 + 9 * u1 * y1m1 - 8 * u3 * (x1my1 + u1 * y2)) + 3 * u2 * (4 * u3m3 * x2 - u1 * y1m1 - 2 * u3 * (x1my1 + u1 * y2))) - 2 * u2m2 * (y1m2 + 3 * u3 * y3) + u2 * (y1m1m1 + 6 * u3 * y1m2 - 12 * u3m3 * y3) + u3 * (-3 * y1m1m1 + 8 * u3 * y1m2 - 6 * u3m3 * y3));
//				double T2_over = (u2 - u3) * (T1mT1 * u1 * u13 + y1m1 + 2 * T1 * (u3mx1 - u1 * y1) - 2 * u3my2) * (T1mT1 * u1 * u123 + y1m1 + 2 * T1 * (u2 * x1 + u3mx1 - u1 * y1) - 2 * (u2 + u3) * y2);
//				// T2_under * T2 + T2_over == 0
//				double T2 = -T2_over / T2_under;
//				double T3 = (u1 * T1 + u2 * T2 - y1) / (-u3);
//				if (0 <= T3 && T3 <= T2 && T2 <= T1) {
//					if (flag_succeed && T1 <= T[0]) {
//						continue;
//						// throw std::invalid_argument("Multiple solutions.");
//						// exit(1); // TODO
//					}
//					flag_succeed = true;
//					T[0] = T1;
//					T[1] = T2;
//					T[2] = T3;
//					//break;
//				}
//			}
//		}
//		if (flag_succeed) {
//			return T;
//		}
//	}
//	else {
//		double c3 = u1 * (u1m1 + 3 * u1 * u2 + 2 * u2m2);
//		double c2 = (-3 * y1) * u1m1 + (6 * u2 * x1 - 3 * u2 * y1) * u1 + 6 * u2m2 * x1;
//		double c1 = (3 * y1m1 - 6 * u2 * y2) * u1 + 12 * x2 * u2m2 - 6 * x1my1 * u2;
//		double c0 = -12 * y3 * u2m2 + 6 * y2 * u2 * y1 - y1m1m1;
//
//		double coeffs[4] = { c0, c1, c2, c3 };
//		static double x[3];
//		int num_root = solve_polynomial_real(6, coeffs, x);
//		bool flag_succeed = false;
//		double* T = new double[3]();
//		for (int i = 0; i < num_root; ++i) {
//			if (x[i] > 0) {
//				double T1 = x[i];
//				double D2 = 2 * u2m2;
//				double D1 = 2 * T1 * u1 * u2 - 2 * u2 * y1;
//				double D0 = T1 * T1 * u1 * u12 - 2 * T1 * (u1 * y1 + u2 * x1) + y1m1 - 2 * u2 * y2;
//				double Delta = D1 * D1 - 4 * D2 * D0;
//				if (Delta < 0) {
//					continue;
//				}
//				double T2 = (-D1 + sqrt(Delta)) / (2 * D2);
//				double T3 = (u1 * T1 + u2 * T2 - y1) / (-u2);
//				if (0 <= T3 && T3 <= T2 && T2 <= T1) {
//					if (flag_succeed) {
//						throw "Multiple solutions.";
//						exit(1); // TODO
//					}
//					flag_succeed = true;
//					T[0] = T1;
//					T[1] = T2;
//					T[2] = T3;
//					//break;
//				}
//				T2 = (-D1 - sqrt(Delta)) / (2 * D2);
//				T3 = (u1 * T1 + u2 * T2 - y1) / (-u2);
//				if (0 <= T3 && T3 <= T2 && T2 <= T1) {
//					if (flag_succeed) {
//						throw "Multiple solutions.";
//						exit(1); // TODO
//					}
//					flag_succeed = true;
//					T[0] = T1;
//					T[1] = T2;
//					T[2] = T3;
//					//break;
//				}
//			}
//		}
//		if (flag_succeed) {
//			return T;
//		}
//	}
//	//throw "No solutions.";
//	//exit(1); // TODO
//	return nullptr;
//}

// Solve tangent markers in 3rd-order problems based on Groebner basis.
// If positive_time is true, then x0->arc1->arc2->(3,2); else, (3,2)->arc2->arc1->x0.
double* Planner_MIM::solution_2arc_tangent_3rd_order(double* x0, const arc& arc1, const arc& arc2, const Constraint& constraint, bool positive_time) {
	double xf3 = (positive_time == arc2.sign) ? constraint.M_max[3] : constraint.M_min[3];
	double u1 = constraint.cal_u(arc1);
	double u2 = constraint.cal_u(arc2) - u1;
	double x1 = x0[0];
	double x2 = x0[1];
	double y2 = -x0[1];
	double y3 = xf3 - x0[2];

	double c6 = u1 * u1 * (u1 + u2);
	double c5 = 6 * u1 * x1 * (u1 + u2);
	double c4 = (-6 * y2) * u1 * u1 + (12 * x1 * x1 + 12 * u2 * x2) * u1 + 9 * u2 * x1 * x1;
	double c3 = (-12 * u2 * y3 - 24 * x1 * y2) * u1 + 8 * x1 * x1 * x1 + 36 * u2 * x2 * x1;
	double c2 = -24 * x1 * x1 * y2 - 36 * u2 * y3 * x1 + 36 * u2 * x2 * x2 + 12 * u1 * y2 * y2;
	double c1 = 24 * x1 * y2 * y2 - 72 * u2 * x2 * y3;
	double c0 = -8 * y2 * y2 * y2 + 36 * u2 * y3 * y3;

	double coeffs[7] = { c0, c1, c2, c3, c4, c5, c6 };
	static double x[6];
	int num_root = solve_polynomial_real(6, coeffs, x);

	bool flag_succeed = false;
	double* T = new double[2]();
	for (int i = 0; i < num_root; ++i) {
		if (positive_time == x[i] > 0) {
			//  T2^2  == (2*y2 - u1*T1^2 - 2*x1*T1)/u2
			double T1 = x[i];
			double T2 = (2.0 * y2 - u1 * T1 * T1 - 2 * x1 * T1) / u2;
			if (T2 < 0) {
				continue;
			}
			T2 = positive_time ? sqrt(T2) : -sqrt(T2);
			double xf1 = x1 + u1 * T1 + u2 * T2;
			if (((positive_time && 0 <= T2 && T2 <= T1) || ((!positive_time) && 0 >= T2 && T2 >= T1)) && (xf3 * xf1 <= 0)) {
				double x1[3];
				dynamics_onestep(3, x0, x1, u1, T1 - T2);
				double x2[3];
				dynamics_onestep(3, x1, x2, u2 + u1, T2);
				if (fabs(x2[2] - xf3) < constraint.EPSILON && fabs(x2[1]) < constraint.EPSILON && (positive_time == arc2.sign ? x2[0] <= 0 : x2[0] >= 0)) {
					if (flag_succeed) {
						throw std::invalid_argument("Multiple solutions.");
						exit(1); // TODO
					}
					flag_succeed = true;
					T[0] = T1;
					T[1] = T2;
					//break;
				}
			}
		}
	}
	if (flag_succeed) {
		return T;
	}
	delete[] T;
	return nullptr;
}

int Planner_MIM::solve_polynomial_real(int order, double const* coeffs, double* x) {
	order++;
	while (fabs(coeffs[--order]) < 1e-6);
	if (order <= 0) {
		return 0;
	}
	else if (order == 1) {
		x[0] = -coeffs[0] / coeffs[1];
		return 1;
	}
	else if (order == 2) {
		double Delta = coeffs[1] * coeffs[1] - 4 * coeffs[0] * coeffs[2];
		if (fabs(Delta) < 1E-6) {
			x[0] = -coeffs[1] / (2 * coeffs[2]);
			x[1] = x[0];
			return 1;
		}
		else if (Delta < 0) {
			return 0;
		}
		else {
			double sqrtDelta = sqrt(Delta);
			x[0] = (-coeffs[1] + sqrtDelta) / (2 * coeffs[2]);
			x[1] = (-coeffs[1] - sqrtDelta) / (2 * coeffs[2]);
			return 2;
		}
	}
	else if (order == 3) {
		// x^3+A*x^2+B*x+C=0
		double A = coeffs[2] / coeffs[3];
		double B = coeffs[1] / coeffs[3];
		double C = coeffs[0] / coeffs[3];
		// Let x = y - A/3, then y^3 + p y + q = 0
		double p = B - A * A / 3;
		double q = 2 * A * A * A / 27 - A * B / 3 + C;
		double discriminant = (q * q) / 4 + (p * p * p) / 27;
		if (discriminant > 0) { // One real root
			double sqrt_D = sqrt(discriminant);
			double u = std::cbrt(-q / 2 + sqrt_D);
			double v = std::cbrt(-q / 2 - sqrt_D);
			double y = u + v;
			x[0] = y - A / 3;
			return 1;
		}
		else if (abs(discriminant) < 1e-6) { // All real roots, at least two are equal
			double u = cbrt(-q / 2);
			x[0] = 2 * u - A / 3;
			x[0] = -u - A / 3;
			return 2;
		}
		else { // Three real roots, all real and distinct
			const static double PI = 3.14159265358979323846;
			double r = sqrt(-p * p * p / 27);
			double phi = acos(-q / (2 * sqrt(-p * p * p / 27)));
			double t = 2 * sqrt(-p / 3);
			x[0] = t * cos(phi / 3) - A / 3;
			x[1] = t * cos((phi + 2 * PI) / 3) - A / 3;
			x[2] = t * cos((phi + 4 * PI) / 3) - A / 3;
			return 3;
		}
	}
	else if (order == 4) { // Use https://github.com/sasamil/Quartic.git
		// x^4 + A*x^3 + B*x^2 + C*x + D = 0
		double a = coeffs[3] / coeffs[4];
		double b = coeffs[2] / coeffs[4];
		double c = coeffs[1] / coeffs[4];
		double d = coeffs[0] / coeffs[4];


		// cubic resolvent
		// y^3 - b*y^2 + (ac-4d)*y - a^2*d - c^2+4*b*d = 0
		double coeffs[4] = { -a * a * d - c * c + 4. * b * d ,a * c - 4. * d,-b,1.0 };
		static double x3[3];
		int iZeroes = solve_polynomial_real(3, coeffs, x3);

		double q1, q2, p1, p2, sqD;
		double y = x3[0];
		if (iZeroes != 1) { // THE ESSENCE - choosing Y with maximal absolute value !
			if (fabs(x3[1]) > fabs(y)) y = x3[1];
			if (fabs(x3[2]) > fabs(y)) y = x3[2];
		}

		// h1+h2 = y && h1*h2 = d  <=>  h^2 -y*h + d = 0    (h === q)
		double D = y * y - 4 * d;
		if (fabs(D) < 1E-6) //in other words - D==0
		{
			q1 = q2 = y * 0.5;
			// g1+g2 = a && g1+g2 = b-y   <=>   g^2 - a*g + b-y = 0    (p === g)
			D = a * a - 4 * (b - y);
			if (fabs(D) < 1E-6) { //in other words - D==0
				p1 = p2 = a * 0.5;
			}
			else {
				sqD = sqrt(D);
				p1 = (a + sqD) * 0.5;
				p2 = (a - sqD) * 0.5;
			}
		}
		else
		{
			sqD = sqrt(D);
			q1 = (y + sqD) * 0.5;
			q2 = (y - sqD) * 0.5;
			// g1+g2 = a && g1*h2 + g2*h1 = c       ( && g === p )  Krammer
			p1 = (a * q1 - c) / (q1 - q2);
			p2 = (c - a * q2) / (q1 - q2);
		}
		int num_root = 0;
		// solving quadratic eq. - x^2 + p1*x + q1 = 0
		D = p1 * p1 - 4 * q1;
		if (fabs(D) < 1E-12) {
			x[num_root++] = -0.5 * p1;
		}
		else if (D > 0.0) {
			sqD = sqrt(D);
			x[num_root++] = (-p1 + sqD) * 0.5;
			if (fabs(sqD) > 1E-6) {
				x[num_root++] = (-p1 - sqD) * 0.5;
			}
		}

		// solving quadratic eq. - x^2 + p2*x + q2 = 0
		D = p2 * p2 - 4 * q2;
		if (fabs(D) < 1E-12) {
			x[num_root++] = -0.5 * p2;
		}
		else if (D > 0.0) {
			sqD = sqrt(D);
			x[num_root++] = (-p2 + sqD) * 0.5;
			if (fabs(sqD) > 1E-6) {
				x[num_root++] = (-p2 - sqD) * 0.5;
			}
		}

		return num_root;
	}
	else {
		static double z[12];
		gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(order + 1);
		gsl_poly_complex_solve(coeffs, order + 1, w, z);
		gsl_poly_complex_workspace_free(w);
		int num = 0;
		for (int i = 0; i < order; ++i) {
			if (fabs(z[2 * i + 1]) < 1e-10) {
				x[num++] = z[2 * i];
			}
		}
		return num;
	}
}

}