#pragma once
#include "arcs.h"
#include <vector>
using namespace std;

class Planner {
public:
	static vector<arc> plan(int order, double* x0, double* xf, Constraint& constraint, bool flag_consider_position);
	static waypoints end_points(int order, const double* x0, const Constraint& constraint, const arc* arcs, int length);
	static waypoints start_points(int order, const double* xf, const Constraint& constraint, const arc* arcs, int length);
	static void dynamics_onestep(int order, const double* x0, double* xf, double u, double dt);
private:
	static inline arc plan_1st_order(double x0, double xf, const Constraint& constraint);
	static vector<arc> plan_2nd_order(double* x0, double* xf, Constraint& constraint, bool flag_consider_position, int direction);
	static vector<arc> plan_3rd_order(double* x0, double* xf, Constraint& constraint, bool flag_consider_position, int direction);
	static bool feasible(int order, const double* x0, const Constraint& constraint, const arc* arcs, int length);
	inline static void change_sign_x(double* x, int order);
	inline static void change_sign_arcs(arc* arcs, int length);
	static double norm_inf(double* x, double* y, int n);
	static double proper_position(int order, double* x0, double* xf, Constraint& constraint);
	static double* solution_3arc_3rd_order(double* x0, double* xf, const arc& arc1, const arc& arc2, const arc& arc3, const Constraint& constraint);
	static double* solution_2arc_tangent_3rd_order(double* x0, const arc& arc1, const arc& arc2, const Constraint& constraint, bool positive_time);
};