#include "demo.h"
#include "Planner.h"
#include "Interpolator.h"
#include <iostream>
#include <chrono>
using namespace std;

void demo_3rd_order_00_3_2_000() {
	int order = 3;
	double M_max[4] = { 1.0, 1.0, 1.5, 4.0 };
	double M_min[4] = { -1.0, -1.0, -1.5, -4.0 };

	Constraint constraint;
	constraint.copy(order, M_max, M_min);
	double x0[3] = { 1,-0.375, 3.999 };
	double xf[3] = { 0, 0, 4 };

	chrono::steady_clock::time_point start, end;
	start = chrono::high_resolution_clock::now();

	vector<arc> arcs = Planner::plan(order, x0, xf, constraint, true);

	end = chrono::high_resolution_clock::now();
	printf("Time: %.3lfms.\n", (end - start).count() * 1.0e-6); // 0.074 ms

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

	double Ts = 1e-3;
	Interpolator interpolator(order);
	interpolator.interpolate(x0, arcs.data(), arcs.size(), constraint, Ts);
	interpolator.write_csv(R"(..\data\3rd_order\00_3_2_000.csv)");

}

void demo_3rd_order_0102010() {
	int order = 3;
	double M_max[4] = { 1.0, 1.0, 1.5, 4.0 };
	double M_min[4] = { -1.0, -1.0, -1.5, -4.0 };

	Constraint constraint;
	constraint.copy(order, M_max, M_min);
	double x0[3] = { 1, -3.0 / 8.0, 4 };
	double xf[3] = { -0.1, 0.1, -1 };

	chrono::steady_clock::time_point start, end;
	start = chrono::high_resolution_clock::now();

	vector<arc> arcs = Planner::plan(order, x0, xf, constraint, true);

	end = chrono::high_resolution_clock::now();
	printf("Time: %.3lfms.\n", (end - start).count() * 1.0e-6); // 0.022 ms

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

	double Ts = 1e-3;
	Interpolator interpolator(order);
	interpolator.interpolate(x0, arcs.data(), arcs.size(), constraint, Ts);
	interpolator.write_csv(R"(..\data\3rd_order\0102010.csv)");
}