#include "demo.h"
#include "Planner.h"
#include <iostream>
#include <chrono>
using namespace std;

void demo_3rd_order() {
	int order = 3;
	double M_max[4] = { 1.0, 1.0, 1.5, 4.0 };
	double M_min[4] = { -1.0, -1.0, -1.5, -4.0 };

	Constraint constraint;
	constraint.copy(order, M_max, M_min);

	// No tangent markers
	//double x0[3] = { 1, -3.0 / 8.0, 4 };
	//double x0[3] = { -1, -1, -5.0 / 6.0 };
	//double x0[3] = { 0,-0.1, -1.9 };
	//double xf[3] = { 0, 0, 0 };

	// tangent marker
	double x0[3] = { 1,-0.375, 3.999 };
	double xf[3] = { 0, 0, 4 };

	chrono::steady_clock::time_point start, end;
	start = chrono::high_resolution_clock::now();

	vector<arc> arcs = Planner::plan(order, x0, xf, constraint, true);

	end = chrono::high_resolution_clock::now();
	printf("Time: %.3lfms.\n", (end - start).count() * 1.0e-6); // ms

	for (int i = 0; i < arcs.size(); i++) {
		cout << "arc " << i << ": order=" << arcs[i].order << ", sign=" << arcs[i].sign << ", tangent=" << arcs[i].tangent << ", time=" << arcs[i].time << endl;
	}
	waypoints xs = Planner::end_points(3, x0, constraint, arcs.data(), arcs.size());
	for (int i = 0; i <= arcs.size(); i++) {
		cout << "x" << i << ": ";
		for (int j = 0; j < order; j++) {
			cout << xs[i][j] << " ";
		}
		cout << endl;
	}
}

void demo_2nd_order() {
	int order = 2;
	double M_max[3] = { 1.0, 1.0, 1.5 };
	double M_min[3] = { -1.0, -1.0, -1.5 };

	Constraint constraint;
	constraint.copy(order, M_max, M_min);

	double x0[2] = { 0, 0.5 };
	double xf[2] = { 0, 1.5 };

	chrono::steady_clock::time_point start, end;
	start = chrono::high_resolution_clock::now();

	vector<arc> arcs = Planner::plan(order, x0, xf, constraint, true);

	end = chrono::high_resolution_clock::now();
	printf("Time: %.3lfms.\n", (end - start).count() * 1.0e-6); // ms

	for (int i = 0; i < arcs.size(); i++) {
		cout << "arc " << i << ": order=" << arcs[i].order << ", sign=" << arcs[i].sign << ", tangent=" << arcs[i].tangent << ", time=" << arcs[i].time << endl;
	}
}

void demo_end_point() {
	int order = 2;
	double M_max[3] = { 1.0, 1.0, 1.5 };
	double M_min[3] = { -1.0, -1.0, -1.5 };

	Constraint constraint;
	constraint.copy(order, M_max, M_min);

	double x0[2] = { 0, 0 };
	vector<arc> arcs;
	arcs.push_back(arc(0, true, 0, 1.0));
	arcs.push_back(arc(1, true, 0, 0.5));
	arcs.push_back(arc(0, false, 0, 1.0));

	waypoints xs = Planner::end_points(order, x0, constraint, arcs.data(), arcs.size());
	for (int i = 0; i <= arcs.size(); i++) {
		cout << "x" << i << ": ";
		for (int j = 0; j < order; j++) {
			cout << xs[i][j] << " ";
		}
		cout << endl;
	}
}