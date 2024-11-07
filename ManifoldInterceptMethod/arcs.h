#pragma once
#include <algorithm>
using namespace std;

struct arc{
	arc();
	arc(int order_, bool sign_, int tangent_, double time_ = 0);

	int order;
	bool sign; // true: +; false: -.
	int tangent; // 0: arc; k>0: tangent marker
	double time; // >=0
};

class Constraint {
public:
	void copy(int order, const double* M_max, const double* M_min);
	void clear();
	~Constraint();
	inline double cal_u(const arc& arc) const;
	void change_sign(int order);
public:
	double* M_max = nullptr;
	double* M_min = nullptr;
	double EPSILON = 1e-6;
};

// Calculate the control input u for the arc.
inline double Constraint::cal_u(const arc& arc) const {
	return arc.order == 0 ? (arc.sign ? M_max[0] : M_min[0]) : 0;
}

struct waypoints {
	waypoints();
	int order;
	int length;
	double* xs; // xs[i*order+j]: xs[i][j]
	void clear();
	void set(int order, int length);
	inline double* operator[](int i);
	~waypoints();
};

inline double* waypoints::operator[](int i) {
	return xs + i * order;
}