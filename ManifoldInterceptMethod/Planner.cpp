#include "Planner.h"
#include <gsl/gsl_poly.h>
#include <stdexcept>
#include <iostream>
using namespace std;

// Plan a 1st-order trajectory from x0 to xf under constraint. (trivial)
inline arc Planner::plan_1st_order(double x0, double xf, const Constraint& constraint) {
	return x0 < xf ? arc(0, true, 0, (xf - x0) / constraint.M_max[0]) : arc(0, false, 0, (xf - x0) / constraint.M_min[0]);
}

// Plan a 2nd-order trajectory from x0 to xf under constraint.
// Only: 00, 010
// If direction is 1/-1, then directly consider the case where x0 is lower/higher than xf.
vector<arc> Planner::plan_2nd_order(double* x0, double* xf, Constraint& constraint, bool flag_consider_position, int direction) {
	if (direction == 0) {
		arc arc_try = plan_1st_order(x0[0], xf[0], constraint);
		waypoints xs_try = end_points(2, x0, constraint, &arc_try, 1);
		if (norm_inf(xs_try[xs_try.length - 1], xf, 2) < constraint.EPSILON) {
			vector<arc> result;
			result.push_back(arc_try);
			if (flag_consider_position && (!feasible(2, x0, constraint, result.data(), result.size()))) {
				result.clear();
			}
			return result;
		}
		direction = xs_try[xs_try.length - 1][1] < xf[1] ? 1 : -1;
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
		result.push_back(arc(0, false, 0, (constraint.M_max[1] - xf[0]) / constraint.M_max[0]));
		if (flag_consider_position && !feasible(2, x0, constraint, result.data(), result.size())) {
			result.clear();
		}
		return result;
	}
	// Consider 00
	// coeffs[0] + coeffs[1]* T2 + coeffs[2] * T2^2 = 0
	double coeffs_00[3] = { xf[0] * xf[0] - x0[0] * x0[0] + 2.0 * constraint.M_max[0] * (x0[1] - xf[1]),2.0 * xf[0] * (constraint.M_max[0] - constraint.M_min[0]),constraint.M_min[0] * (constraint.M_min[0] - constraint.M_max[0]) };
	gsl_poly_complex_workspace* w_00 = gsl_poly_complex_workspace_alloc(3);
	double z_00[4];
	gsl_poly_complex_solve(coeffs_00, 3, w_00, z_00);
	gsl_poly_complex_workspace_free(w_00);
	double T2 = max(z_00[0], z_00[2]); // T2 > 0 holds
	double T1 = (xf[0] - x0[0] - (constraint.M_min[0] - constraint.M_max[0]) * T2) / constraint.M_max[0]; // T1 > T2 holds
	result.push_back(arc(0, true, 0, T1 - T2));
	result.push_back(arc(0, false, 0, T2));
	if (flag_consider_position && !feasible(2, x0, constraint, result.data(), result.size())) {
		result.clear();
	}
	return result;
}

// Change the sign of x.
inline void Planner::change_sign_x(double* x, int order) {
	for (int k = 0; k < order; k++) {
		x[k] = -x[k];
	}
}

inline void Planner::change_sign_arcs(arc* arcs, int length) {
	for (int i = 0; i < length; ++i) {
		arcs[i].sign = !arcs[i].sign;
	}
}

vector<arc> Planner::plan(int order, double* x0, double* xf, Constraint& constraint, bool flag_consider_position) {
	if (order == 1) {
		vector<arc> result;
		result.push_back(plan_1st_order(x0[0], xf[0], constraint));
		return result;
	}
	if (order == 2) {
		return plan_2nd_order(x0, xf, constraint, flag_consider_position, 0);
	}
	return plan_3rd_order(x0, xf, constraint, flag_consider_position, 0);
}

// Calculate the state vector at the end of the arc.
void Planner::dynamics_onestep(int order, const double* x0, double* xf, double u, double dt) {
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

// Calculate the end points of the trajectory from the initial state vector x0.
waypoints Planner::end_points(int order, const double* x0, const Constraint& constraint, const arc* arcs, int length) {
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
	//double* T = new double[order + 1]; // T[k] = t^k/k!
	for (int i = 0, i_arc = 0; i < N_real; i++, i_arc++) {
		while (arcs[i_arc].tangent != 0) {
			i_arc++;
		}
		// x moves from x[k] to x[k+1] along arcs[i_arc]
		dynamics_onestep(order, result[i], result[i + 1], constraint.cal_u(arcs[i_arc]), arcs[i_arc].time);
	}
	//delete[] T;
	return result;
}

// Calculate the end points of the trajectory to the terminal state vector xf.
waypoints Planner::start_points(int order, const double* xf, const Constraint& constraint, const arc* arcs, int length) {
	arc* arcs_reverse = new arc[length];
	for (int i = 0; i < length; i++) {
		arcs_reverse[i] = arcs[length - i - 1];
		arcs_reverse[i].time = -arcs_reverse[i].time;
	}
	waypoints result = end_points(order, xf, constraint, arcs_reverse, length);
	delete[] arcs_reverse;
	return result;
}


// Check if the trajectory is feasible.
bool Planner::feasible(int order, const double* x0, const Constraint& constraint, const arc* arcs, int length) {
	double* x = new double[order];
	for (int k = 0; k < order; k++) {
		if ((!isinf(constraint.M_min[k + 1]) && x0[k] + constraint.EPSILON < constraint.M_min[k + 1]) || (!isinf(constraint.M_max[k + 1]) && x0[k] - constraint.EPSILON > constraint.M_max[k + 1])) {
			delete[] x;
			return false;
		}
		x[k] = x0[k];
	}
	double* x_old = new double[order];
	double* coeffs = new double[order]();
	double* T = new double[order](); // T[k] = t^k/k!
	double* z = new double[order + 1]();
	for (int i = 0; i < length; i++) {
		if (arcs[i].tangent != 0) {
			continue;
		}
		double u = constraint.cal_u(arcs[i]);
		// Check the local minimum or maximum of x[k] along the trajectory.
		T[0] = 1;
		for (int k = 1; k < order; k++) {
			T[k] = T[k - 1] * arcs[i].time / k;
			if (isinf(constraint.M_max[k]) && isinf(constraint.M_min[k])) {
				continue;
			}
			// Find the root of x[k-1]: $k$th-order polynomial
			for (int j = 0; j < k; j++) {
				coeffs[j] = T[j] * x[k - j - 1];
			}
			coeffs[k] = T[k] * u;
			int k1 = k;
			while (fabs(coeffs[k1]) < 1e-10 && k1 >= 0) {
				k1--;
			}
			if (k1 <= 0) {
				continue;
			}
			gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(k1 + 1);
			gsl_poly_complex_solve(coeffs, k1 + 1, w, z);
			gsl_poly_complex_workspace_free(w);
			for (int l = 0; l < k1; l++) {
				if (fabs(z[2 * l + 1]) < constraint.EPSILON && z[2 * l] > 0 && z[2 * l] < arcs[i].time) {
					// A local minimum or maximum of x[k] is found.
					double xk = x[k];
					double Tk = 1;
					for (int j = 1; j <= k; j++) {
						Tk *= z[2 * l] / j;
						xk += Tk * x[k - j];
					}
					if ((!isinf(constraint.M_min[k + 1]) && xk + constraint.EPSILON < constraint.M_min[k + 1]) || (!isinf(constraint.M_max[k + 1]) && xk - constraint.EPSILON > constraint.M_max[k + 1])) {
						delete[] x;
						delete[] x_old;
						delete[] coeffs;
						delete[] z;
						return false;
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
			if ((!isinf(constraint.M_min[k + 1]) && x[k] + constraint.EPSILON < constraint.M_min[k + 1]) || (!isinf(constraint.M_max[k + 1]) && x[k] - constraint.EPSILON > constraint.M_max[k + 1])) {
				delete[] x;
				delete[] x_old;
				delete[] coeffs;
				delete[] z;
				return false;
			}
		}
	}
	delete[] x;
	delete[] x_old;
	delete[] coeffs;
	delete[] z;
	return true;
}


// Plan a 3rd-order trajectory from x0 to xf under constraint.
// First Step: 00, 010
// Second Step: s00, s010
// If direction is 1/-1, then directly consider the case where x0 is lower/higher than xf.
vector<arc> Planner::plan_3rd_order(double* x0, double* xf, Constraint& constraint, bool flag_consider_position, int direction) {
	if (direction == 0) {
		vector<arc> arc_try = plan_2nd_order(x0, xf, constraint, false, 0);
		if (arc_try.empty()) {
			return arc_try;
		}
		waypoints xs_try = end_points(3, x0, constraint, arc_try.data(), arc_try.size());
		if (norm_inf(xs_try[xs_try.length - 1], xf, 3) < constraint.EPSILON) {
			return arc_try;
		}
		direction = xs_try[xs_try.length - 1][2] < xf[2] ? 1 : -1;
	}
	if (direction < 0) { // x0 is higher than xf. It must be feasible.
		constraint.change_sign(3);
		change_sign_x(x0, 3);
		change_sign_x(xf, 3);
		vector<arc> result = plan_3rd_order(x0, xf, constraint, flag_consider_position, 1);
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
		double propos2 = proper_position(3, maxV, xf, constraint);
		waypoints xs2 = start_points(3, xf, constraint, result2.data(), result2.size());
		//printf("propos2=%lf, xs_before=%lf, xs2=%lf\n", propos2, xs1[xs1.length - 1][2], xs2[xs2.length - 1][2]);
		if (xs1[xs1.length - 1][2] <= xf[2] - propos2) {
			result.push_back(arc(2, true, 0, (xf[2] - propos2 - xs1[xs1.length - 1][2]) / constraint.M_max[2]));
			for (int i = 0; i < result2.size(); i++) {
				result.push_back(result2[i]);
			}
			xs2.clear();
			xs1.clear();
			return result;
		}
		//for (int i = 0; i < xs1.length; i++) {
		//	cout << "x" << i << ": ";
		//	for (int j = 0; j < xs1.order; j++) {
		//		cout << xs1[i][j] << " ";
		//	}
		//	cout << endl;
		//}
		//printf("propos=%lf, xs_before=%lf, xf=%lf\n", proper_position(3, xs1[0], xf, constraint), xs1[0][2], xf[2]);
		//printf("propos=%lf, xs_before=%lf, xf=%lf\n", proper_position(3, xs1[1], xf, constraint), xs1[1][2], xf[2]);
		//printf("propos=%lf, xs_before=%lf, xf=%lf\n", proper_position(3, xs1[2], xf, constraint), xs1[2][2], xf[2]);
		//printf("propos=%lf, xs_before=%lf, xf=%lf\n", proper_position(3, xs1[3], xf, constraint), xs1[3][2], xf[2]);
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
	// Step 2: Find where the trajectory enters the manifold.
	bool flag_succeed_noposition = false;
	while (result.size() > 1) {
		double propos = proper_position(3, xs_before[xs_before.length - 1], xf, constraint);
		//printf("propos=%lf, xs_before=%lf, xf=%lf\n", propos, xs_before[xs_before.length - 1][2], xf[2]);
		if (xs_before[xs_before.length - 1][2] <= xf[2] - propos) {
			break;
		}
		result.pop_back();
		result[result.size() - 1].time = 0;
		--xs_before.length;
	}
	// xs1[xs1.length-1] enters the manifold along the arc result[result.size()-1].
	// Step 3: Create the whole trajectory. Consider s00, s010
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
		if (feasible(2, xs_before[xs_before.length - 1], constraint, result.data() + (result.size() - 3), 3)) {
			flag_succeed_noposition = true;
		}
		else {
			result.pop_back();
			result.pop_back();
		}
	}
	else {
		throw invalid_argument("No solution for s00.");
		exit(1);
	}
	if (!flag_succeed_noposition) {
		if (isinf(constraint.M_min[1])) {
			throw invalid_argument("No solution for s01.");
			exit(1);
		}
		// Step 3.2: s010
		arc arc4 = plan_1st_order(constraint.M_min[1], xf[0], constraint);
		waypoints xsf = start_points(3, xf, constraint, &arc4, 1);
		arc3 = arc(1, false, 0);
		double* T_s01 = solution_3arc_3rd_order(xs_before[xs_before.length - 1], xsf[xsf.length - 1], result[result.size() - 1], arc2, arc3, constraint);
		result[result.size() - 1].time = T_s01[0] - T_s01[1];
		arc2.time = T_s01[1] - T_s01[2];
		arc3.time = T_s01[2];
		delete[] T_s01;
		result.push_back(arc2);
		result.push_back(arc3);
		if (feasible(2, xs_before[xs_before.length - 1], constraint, result.data() + (result.size() - 3), 3)) {
			flag_succeed_noposition = true;
		}
		else {
			result.pop_back();
			result.pop_back();
		}
	}
	if (!flag_succeed_noposition) {
		throw invalid_argument("No solution.");
		exit(1);
	}
	if (feasible(3, x0, constraint, result.data(), result.size())) {
		return result;
	}
	// Step 4: Consider tangent markers.
	// Only -3 can support tangent markers, i.e., (-3,2).
	result.clear();
	// Step 4.1: like 00(-3,2)000
	// Step 4.1.1: 00(-3,2)*
	arc arc1(0, true, 0);
	arc2 = arc(0, false, 0);
	double* T_tangent = solution_2arc_tangent_3rd_order(x0, arc1, arc2, constraint, true);
	if (T_tangent != nullptr) {
		arc1.time = T_tangent[0] - T_tangent[1];
		arc2.time = T_tangent[1];
		delete[] T_tangent;
		result.push_back(arc1);
		result.push_back(arc2);
		waypoints xs1 = end_points(3, x0, constraint, result.data(), result.size());
		result.push_back(arc(3, true, 2));
		vector<arc> result2 = plan_3rd_order(xs1[2], xf, constraint, false, -1);
		for (int i = 0; i < result2.size(); ++i) {
			result.push_back(result2[i]);
		}
		return result;
	}
	// Step 4.1.2: 010(-3,2)*
	arc1 = plan_1st_order(x0[0], constraint.M_max[1], constraint);
	arc2 = arc(1, true, 0);
	arc3 = arc(0, false, 0);
	waypoints x0now_pre = end_points(3, x0, constraint, &arc1, 1);
	T_tangent = solution_2arc_tangent_3rd_order(x0now_pre[1], arc2, arc3, constraint, true);
	if (T_tangent != nullptr) {
		arc2.time = T_tangent[0] - T_tangent[1];
		arc3.time = T_tangent[1];
		delete[] T_tangent;
		result.push_back(arc1);
		result.push_back(arc2);
		result.push_back(arc3);
		waypoints xs1 = end_points(3, x0, constraint, result.data(), result.size());
		result.push_back(arc(3, true, 2));
		vector<arc> result2 = plan_3rd_order(xs1[2], xf, constraint, false, -1);
		for (int i = 0; i < result2.size(); ++i) {
			result.push_back(result2[i]);
		}
		return result;
	}
	// Step 4.2: like 000(-3,2)00
	// Step 4.2.1: *(-3,2)00
	arc1 = arc(0, false, 0);
	arc2 = arc(0, true, 0);
	T_tangent = solution_2arc_tangent_3rd_order(xf, arc1, arc2, constraint, false);
	if (T_tangent != nullptr) {
		arc1.time = T_tangent[1] - T_tangent[0];
		arc2.time = -T_tangent[1];
		delete[] T_tangent;
		arc arcs2[] = { arc2, arc1 };
		waypoints xs1 = start_points(3, xf, constraint, arcs2, 2);
		vector<arc> result2 = plan_3rd_order(x0, xs1[2], constraint, false, -1);
		for (int i = 0; i < result2.size(); ++i) {
			result.push_back(result2[i]);
		}
		result.push_back(arc(3, true, 2));
		result.push_back(arc2);
		result.push_back(arc1);
		return result;
	}
	// Step 4.2.1: *(-3,2)010
	arc1 = plan_1st_order(constraint.M_max[1], xf[0], constraint);
	arc2 = arc(1, true, 0);
	arc3 = arc(0, true, 0);
	waypoints x0now_post = start_points(3, xf, constraint, &arc1, 1);
	T_tangent = solution_2arc_tangent_3rd_order(x0now_pre[1], arc2, arc3, constraint, false);
	if (T_tangent != nullptr) {
		arc2.time = T_tangent[1] - T_tangent[0];
		arc3.time = -T_tangent[1];
		delete[] T_tangent;
		arc arcs2[] = { arc3, arc2, arc1 };
		waypoints xs1 = start_points(3, xf, constraint, arcs2, 3);
		vector<arc> result2 = plan_3rd_order(x0, xs1[3], constraint, false, -1);
		for (int i = 0; i < result2.size(); ++i) {
			result.push_back(result2[i]);
		}
		result.push_back(arc(3, true, 2));
		result.push_back(arc3);
		result.push_back(arc2);
		result.push_back(arc1);
		return result;
	}
	result.clear();
	return result;
}

// Calculate the infinity norm of x-y.
double Planner::norm_inf(double* x, double* y, int n) {
	double dmax = 0;
	for (int i = 0; i < n; i++) {
		dmax = max(dmax, fabs(x[i] - y[i]));
	}
	return dmax;
}

// Calculate the proper position of the trajectory.
double Planner::proper_position(int order, double* x0, double* xf, Constraint& constraint) {
	vector<arc> arcs = plan(order - 1, x0, xf, constraint, false);
	double* x0_ = new double[order];
	for (int k = 0; k < order - 1; k++) {
		x0_[k] = x0[k];
	}
	x0_[order - 1] = 0;
	waypoints xs = end_points(order, x0_, constraint, arcs.data(), arcs.size());
	delete[] x0_;
	return xs[xs.length - 1][order - 1];
}

// Solve 3arcs in 3rd-order problems based on Groebner basis.
double* Planner::solution_3arc_3rd_order(double* x0, double* xf, const arc& arc1, const arc& arc2, const arc& arc3, const Constraint& constraint) {
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

	if (fabs(u2 - u3) > 1e-6) {
		double c6 = u1 * u1 * (u1 + u2) * (u1 + u3) * (u1 + u2 + u3) * (u1 + u2 + u3);
		double c5 = 6 * u1 * (u1 + u2) * (u1 + u3) * (u2 * x1 + u3 * x1 - u1 * y1) * (u1 + u2 + u3);
		double c4 = (15 * y1 * y1 - 6 * u2 * y2 - 6 * u3 * y2) * u1 * u1 * u1 * u1 + (18 * u2 * y1 * y1 + 18 * u3 * y1 * y1 - 12 * u2 * u2 * y2 - 12 * u3 * u3 * y2 + 24 * u2 * u3 * x2 - 12 * u2 * u3 * y2 - 24 * u2 * x1 * y1 - 24 * u3 * x1 * y1) * u1 * u1 * u1 + (12 * u2 * u2 * x1 * x1 + 12 * u3 * u3 * x1 * x1 + 3 * u2 * u2 * y1 * y1 + 3 * u3 * u3 * y1 * y1 - 6 * u2 * u2 * u2 * y2 - 6 * u3 * u3 * u3 * y2 + 12 * u2 * u3 * x1 * x1 + 36 * u2 * u3 * u3 * x2 + 36 * u2 * u2 * u3 * x2 + 15 * u2 * u3 * y1 * y1 - 6 * u2 * u3 * u3 * y2 - 6 * u2 * u2 * u3 * y2 - 24 * u2 * u2 * x1 * y1 - 24 * u3 * u3 * x1 * y1 - 60 * u2 * u3 * x1 * y1) * u1 * u1 + (12 * x2 * u2 * u2 * u2 * u3 + 12 * u2 * u2 * u2 * x1 * x1 + 24 * x2 * u2 * u2 * u3 * u3 + 24 * u2 * u2 * u3 * x1 * x1 - 30 * y1 * u2 * u2 * u3 * x1 + 12 * x2 * u2 * u3 * u3 * u3 + 24 * u2 * u3 * u3 * x1 * x1 - 30 * y1 * u2 * u3 * u3 * x1 + 12 * u3 * u3 * u3 * x1 * x1) * u1 + 9 * u2 * u2 * u2 * u3 * x1 * x1 + 18 * u2 * u2 * u3 * u3 * x1 * x1 + 9 * u2 * u3 * u3 * u3 * x1 * x1;
		double c3 = (24 * u2 * y1 * y2 - 24 * u2 * u3 * y3 - 20 * y1 * y1 * y1 + 24 * u3 * y1 * y2) * u1 * u1 * u1 + (36 * u2 * x1 * y1 * y1 - 12 * u3 * y1 * y1 * y1 - 36 * u2 * u3 * u3 * y3 - 36 * u2 * u2 * u3 * y3 - 12 * u2 * y1 * y1 * y1 + 36 * u3 * x1 * y1 * y1 - 24 * u2 * u2 * x1 * y2 - 24 * u3 * u3 * x1 * y2 + 24 * u2 * u2 * y1 * y2 + 24 * u3 * u3 * y1 * y2 - 24 * u2 * u3 * x1 * y2 - 72 * u2 * u3 * x2 * y1 + 24 * u2 * u3 * y1 * y2) * u1 * u1 + (12 * u2 * u2 * x1 * y1 * y1 - 12 * u2 * u3 * u3 * u3 * y3 - 12 * u2 * u2 * u2 * u3 * y3 - 24 * u2 * u2 * u2 * x1 * y2 - 24 * u3 * u3 * u3 * x1 * y2 - 24 * u2 * u2 * u3 * u3 * y3 - 4 * u2 * u3 * y1 * y1 * y1 - 24 * u2 * u2 * x1 * x1 * y1 + 12 * u3 * u3 * x1 * y1 * y1 - 24 * u3 * u3 * x1 * x1 * y1 + 72 * u2 * u3 * u3 * x1 * x2 + 72 * u2 * u2 * u3 * x1 * x2 + 48 * u2 * u3 * x1 * y1 * y1 - 24 * u2 * u3 * x1 * x1 * y1 - 12 * u2 * u3 * u3 * x1 * y2 - 36 * u2 * u3 * u3 * x2 * y1 - 12 * u2 * u2 * u3 * x1 * y2 - 36 * u2 * u2 * u3 * x2 * y1 + 12 * u2 * u3 * u3 * y1 * y2 + 12 * u2 * u2 * u3 * y1 * y2) * u1 + 36 * x2 * u2 * u2 * u2 * u3 * x1 + 8 * u2 * u2 * u2 * x1 * x1 * x1 + 72 * x2 * u2 * u2 * u3 * u3 * x1 - 8 * u2 * u2 * u3 * x1 * x1 * x1 - 36 * y1 * u2 * u2 * u3 * x1 * x1 + 36 * x2 * u2 * u3 * u3 * u3 * x1 - 8 * u2 * u3 * u3 * x1 * x1 * x1 - 36 * y1 * u2 * u3 * u3 * x1 * x1 + 8 * u3 * u3 * u3 * x1 * x1 * x1;
		double c2 = (12 * u2 * u2 * y2 * y2 + 72 * y3 * u2 * u3 * y1 + 12 * u2 * u3 * y2 * y2 - 36 * u2 * y1 * y1 * y2 + 12 * u3 * u3 * y2 * y2 - 36 * u3 * y1 * y1 * y2 + 15 * y1 * y1 * y1 * y1) * u1 * u1 + (12 * u2 * u2 * u2 * y2 * y2 + 36 * y3 * u2 * u2 * u3 * y1 - 12 * u2 * u2 * u3 * y2 * y2 - 72 * x2 * u2 * u2 * u3 * y2 - 72 * x1 * y3 * u2 * u2 * u3 - 12 * u2 * u2 * y1 * y1 * y2 + 48 * x1 * u2 * u2 * y1 * y2 + 36 * y3 * u2 * u3 * u3 * y1 - 12 * u2 * u3 * u3 * y2 * y2 - 72 * x2 * u2 * u3 * u3 * y2 - 72 * x1 * y3 * u2 * u3 * u3 - 12 * u2 * u3 * y1 * y1 * y2 + 72 * x2 * u2 * u3 * y1 * y1 + 48 * x1 * u2 * u3 * y1 * y2 + 3 * u2 * y1 * y1 * y1 * y1 - 24 * x1 * u2 * y1 * y1 * y1 + 12 * u3 * u3 * u3 * y2 * y2 - 12 * u3 * u3 * y1 * y1 * y2 + 48 * x1 * u3 * u3 * y1 * y2 + 3 * u3 * y1 * y1 * y1 * y1 - 24 * x1 * u3 * y1 * y1 * y1) * u1 - 36 * y3 * u2 * u2 * u2 * u3 * x1 + 36 * u2 * u2 * u2 * u3 * x2 * x2 - 24 * y2 * u2 * u2 * u2 * x1 * x1 - 72 * y3 * u2 * u2 * u3 * u3 * x1 + 72 * u2 * u2 * u3 * u3 * x2 * x2 + 24 * y2 * u2 * u2 * u3 * x1 * x1 - 72 * u2 * u2 * u3 * x1 * x2 * y1 + 36 * y2 * u2 * u2 * u3 * x1 * y1 + 12 * u2 * u2 * x1 * x1 * y1 * y1 - 36 * y3 * u2 * u3 * u3 * u3 * x1 + 36 * u2 * u3 * u3 * u3 * x2 * x2 + 24 * y2 * u2 * u3 * u3 * x1 * x1 - 72 * u2 * u3 * u3 * x1 * x2 * y1 + 36 * y2 * u2 * u3 * u3 * x1 * y1 + 12 * u2 * u3 * x1 * x1 * y1 * y1 - 12 * u2 * u3 * x1 * y1 * y1 * y1 - 24 * y2 * u3 * u3 * u3 * x1 * x1 + 12 * u3 * u3 * x1 * x1 * y1 * y1;
		double c1 = (72 * y3 * u2 * u2 * u3 * y2 - 24 * u2 * u2 * y1 * y2 * y2 + 72 * y3 * u2 * u3 * u3 * y2 - 72 * y3 * u2 * u3 * y1 * y1 - 24 * u2 * u3 * y1 * y2 * y2 + 24 * u2 * y1 * y1 * y1 * y2 - 24 * u3 * u3 * y1 * y2 * y2 + 24 * u3 * y1 * y1 * y1 * y2 - 6 * y1 * y1 * y1 * y1 * y1) * u1 - 72 * x2 * y3 * u2 * u2 * u2 * u3 + 24 * x1 * u2 * u2 * u2 * y2 * y2 - 144 * x2 * y3 * u2 * u2 * u3 * u3 + 72 * x2 * u2 * u2 * u3 * y1 * y2 + 72 * x1 * y3 * u2 * u2 * u3 * y1 - 24 * x1 * u2 * u2 * u3 * y2 * y2 - 24 * x1 * u2 * u2 * y1 * y1 * y2 - 72 * x2 * y3 * u2 * u3 * u3 * u3 + 72 * x2 * u2 * u3 * u3 * y1 * y2 + 72 * x1 * y3 * u2 * u3 * u3 * y1 - 24 * x1 * u2 * u3 * u3 * y2 * y2 - 24 * x2 * u2 * u3 * y1 * y1 * y1 - 24 * x1 * u2 * u3 * y1 * y1 * y2 + 6 * x1 * u2 * y1 * y1 * y1 * y1 + 24 * x1 * u3 * u3 * u3 * y2 * y2 - 24 * x1 * u3 * u3 * y1 * y1 * y2 + 6 * x1 * u3 * y1 * y1 * y1 * y1;
		double c0 = 36 * u2 * u2 * u2 * u3 * y3 * y3 - 8 * u2 * u2 * u2 * y2 * y2 * y2 + 72 * u2 * u2 * u3 * u3 * y3 * y3 - 72 * u2 * u2 * u3 * y1 * y2 * y3 + 8 * u2 * u2 * u3 * y2 * y2 * y2 + 12 * u2 * u2 * y1 * y1 * y2 * y2 + 36 * u2 * u3 * u3 * u3 * y3 * y3 - 72 * u2 * u3 * u3 * y1 * y2 * y3 + 8 * u2 * u3 * u3 * y2 * y2 * y2 + 24 * u2 * u3 * y1 * y1 * y1 * y3 + 12 * u2 * u3 * y1 * y1 * y2 * y2 - 6 * u2 * y1 * y1 * y1 * y1 * y2 - 8 * u3 * u3 * u3 * y2 * y2 * y2 + 12 * u3 * u3 * y1 * y1 * y2 * y2 - 6 * u3 * y1 * y1 * y1 * y1 * y2 + y1 * y1 * y1 * y1 * y1 * y1;

		double coeffs[7] = { c0, c1, c2, c3, c4, c5, c6 };
		double z[12];
		int n = 7;
		while (fabs(coeffs[--n]) < 1e-6);
		gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(n + 1);
		gsl_poly_complex_solve(coeffs, n + 1, w, z);
		gsl_poly_complex_workspace_free(w);
		bool flag_succeed = false;
		double* T = new double[3]();
		for (int i = 0; i < n; ++i) {
			if (fabs(z[2 * i + 1]) < 1e-10 && z[2 * i] > 0) {
				double T1 = z[2 * i];
				double T2_under = -u2 * (-T1 * T1 * T1 * u1 * u1 * u1 * u2 + 3 * T1 * T1 * T1 * u1 * u1 * u1 * u3 - T1 * T1 * T1 * u1 * u1 * u2 * u2 + 3 * T1 * T1 * T1 * u1 * u1 * u2 * u3 + 4 * T1 * T1 * T1 * u1 * u1 * u3 * u3 + T1 * T1 * T1 * u1 * u2 * u2 * u3 + 2 * T1 * T1 * T1 * u1 * u2 * u3 * u3 + T1 * T1 * T1 * u1 * u3 * u3 * u3 + 3 * T1 * T1 * u1 * u1 * u2 * y1 - 9 * T1 * T1 * u1 * u1 * u3 * y1 + T1 * T1 * u1 * u2 * u2 * y1 - 2 * x1 * T1 * T1 * u1 * u2 * u2 - 3 * T1 * T1 * u1 * u2 * u3 * y1 + 6 * x1 * T1 * T1 * u1 * u2 * u3 - 4 * T1 * T1 * u1 * u3 * u3 * y1 + 8 * x1 * T1 * T1 * u1 * u3 * u3 + 3 * x1 * T1 * T1 * u2 * u2 * u3 + 6 * x1 * T1 * T1 * u2 * u3 * u3 + 3 * x1 * T1 * T1 * u3 * u3 * u3 + 2 * y2 * T1 * u1 * u2 * u2 - 6 * y2 * T1 * u1 * u2 * u3 - 3 * T1 * u1 * u2 * y1 * y1 - 8 * y2 * T1 * u1 * u3 * u3 + 9 * T1 * u1 * u3 * y1 * y1 + 6 * x2 * T1 * u2 * u2 * u3 + 2 * x1 * T1 * u2 * u2 * y1 + 12 * x2 * T1 * u2 * u3 * u3 - 6 * x1 * T1 * u2 * u3 * y1 + 6 * x2 * T1 * u3 * u3 * u3 - 8 * x1 * T1 * u3 * u3 * y1 - 6 * y3 * u2 * u2 * u3 - 2 * y2 * u2 * u2 * y1 - 12 * y3 * u2 * u3 * u3 + 6 * y2 * u2 * u3 * y1 + u2 * y1 * y1 * y1 - 6 * y3 * u3 * u3 * u3 + 8 * y2 * u3 * u3 * y1 - 3 * u3 * y1 * y1 * y1);
				double T2_over = (u2 - u3) * (T1 * T1 * u1 * u1 + u3 * T1 * T1 * u1 - 2 * T1 * u1 * y1 + 2 * u3 * x1 * T1 + y1 * y1 - 2 * u3 * y2) * (y1 * y1 - 2 * u3 * y2 - 2 * u2 * y2 + T1 * T1 * u1 * u1 + 2 * T1 * u2 * x1 + 2 * T1 * u3 * x1 - 2 * T1 * u1 * y1 + T1 * T1 * u1 * u2 + T1 * T1 * u1 * u3);
				// T2_under * T2 + T2_over == 0
				double T2 = -T2_over / T2_under;
				double T3 = (u1 * T1 + u2 * T2 - y1) / (-u3);
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
	else {
		double c3 = u1 * (u1 * u1 + 3 * u1 * u2 + 2 * u2 * u2);
		double c2 = (-3 * y1) * u1 * u1 + (6 * u2 * x1 - 3 * u2 * y1) * u1 + 6 * u2 * u2 * x1;
		double c1 = (3 * y1 * y1 - 6 * u2 * y2) * u1 + 12 * x2 * u2 * u2 - 6 * x1 * y1 * u2;
		double c0 = -12 * y3 * u2 * u2 + 6 * y2 * u2 * y1 - y1 * y1 * y1;

		double coeffs[4] = { c0, c1, c2, c3 };
		double z[6];
		int n = 4;
		while (fabs(coeffs[--n]) < 1e-6);
		gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(n + 1);
		gsl_poly_complex_solve(coeffs, n + 1, w, z);
		gsl_poly_complex_workspace_free(w);
		bool flag_succeed = false;
		double* T = new double[3]();
		for (int i = 0; i < n; ++i) {
			if (fabs(z[2 * i + 1]) < 1e-10 && z[2 * i] > 0) {
				double T1 = z[2 * i];
				double D2 = 2 * u2 * u2;
				double D1 = 2 * T1 * u1 * u2 - 2 * u2 * y1;
				double D0 = T1 * T1 * u1 * u1 + u2 * T1 * T1 * u1 - 2 * T1 * u1 * y1 + 2 * u2 * x1 * T1 + y1 * y1 - 2 * u2 * y2;
				double Delta = D1 * D1 - 4 * D2 * D0;
				if (Delta < 0) {
					continue;
				}
				double T2 = (-D1 + sqrt(Delta)) / (2 * D2);
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
	throw "No solutions.";
	exit(1); // TODO
	return nullptr;
}

// Solve tangent markers in 3rd-order problems based on Groebner basis.
// If positive_time is true, then x0->arc1->arc2->(3,2); else, (3,2)->arc2->arc1->x0.
double* Planner::solution_2arc_tangent_3rd_order(double* x0, const arc& arc1, const arc& arc2, const Constraint& constraint, bool positive_time) {
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
	double z[12];
	int n = 7;
	while (fabs(coeffs[--n]) < 1e-6);
	gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(n + 1);
	gsl_poly_complex_solve(coeffs, n + 1, w, z);
	gsl_poly_complex_workspace_free(w);
	bool flag_succeed = false;
	double* T = new double[2]();
	for (int i = 0; i < n; ++i) {
		if (fabs(z[2 * i + 1]) < 1e-10 && (positive_time == z[2 * i] > 0)) {
			//  T2^2  == (2*y2 - u1*T1^2 - 2*x1*T1)/u2
			double T1 = z[2 * i];
			double T2 = (2.0 * y2 - u1 * T1 * T1 - 2 * x1 * T1) / u2;
			if (T2 < 0) {
				continue;
			}
			T2 = positive_time ? sqrt(T2) : -sqrt(T2);
			double xf1 = x1 + u1 * T1 + u2 * T2;
			if (((positive_time && 0 <= T2 && T2 <= T1) || ((!positive_time) && 0 >= T2 && T2 >= T1)) && (xf3 * xf1 <= 0)) {
				if (flag_succeed) {
					throw "Multiple solutions.";
					exit(1); // TODO
				}
				flag_succeed = true;
				T[0] = T1;
				T[1] = T2;
				//break;
			}
		}
	}
	if (flag_succeed) {
		return T;
	}
	delete[] T;
	return nullptr;
}