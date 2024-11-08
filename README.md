# Manifold-Intercept Method (MIM): Plan Time-Optimal Trajectory under $n$th-Order Box-Constraints

[![License](https://img.shields.io/github/license/WangY18/ManifoldInterceptMethod)](https://www.gnu.org/licenses/gpl-3.0.en.html)

<p align="center">
	<img src="https://github.com/user-attachments/assets/62dd3403-aa82-42ab-b190-e46eaf77a6b6" alt="first_show" width="500" />
</p>
<p align="center">
	<b>Figure 1.</b> A fifth-order quasi-optimal trajectory.
</p>

$$
\text{The augmented switching law: }\underline{0}\overline{0}\left(\overline{3},2\right)\overline{0}\underline{0}\overline{03}\underline{01}\overline{0}\underline{2}\overline{01}\underline{0}\overline{4}\underline{01}\overline{0}\underline{2}\overline{01}\underline{0}\underline{3}\overline{01}\underline{0}\overline{0}
$$

## Please Cite:

```
Wang, Yunan, Hu, Chuxiong, Li, Zeyang, Lin, Shize, He, Suqin, & Zhu, Yu (2024). Time-optimal control for high-order chain-of-integrators systems with full state constraints and arbitrary terminal states. IEEE Transactions on Automatic Control.
```

```
@article{wang2024time,
  title={Time-Optimal Control for High-Order Chain-of-Integrators Systems with Full State Constraints and Arbitrary Terminal States},
  author={Wang, Yunan and Hu, Chuxiong and Li, Zeyang and Lin, Shize and He, Suqin and Zhu, Yu},
  journal={IEEE Transactions on Automatic Control},
  year={2024}
}
```

## What is MIM?

MIM is proposed for planning high-order time-optimal trajectories. For example, in Figure 1, we can plan a fifth-order quasi-optimal trajectory, i.e., the position, velocity, acceleration, jerk, snap, and crackle are all bounded. 

MIM can plan:

+ **Time-optimal** trajectories of order $n\leq3$. In other words, the position, velocity, acceleration, and jerk is bounded. Note that only ruckig pro can achieve full constraints which is not open-source. **Other methods all fail to deal with position constraints.**
+ **Quasi-optimal** trajectories of order $n\geq4$.

The advantages of MIM are as follows: 

+ **Asymmetric constraints**. For example, we can require that $-1\leq x_3\leq +\infty$, $-2\leq x_2\leq 3$, $-\infty\leq x_1\leq +\infty$, $-1.5\leq u\leq 1$ in a 3rd-order problem, where $x_3$ is position, $x_2$ is velocity, $x_1$ is acceleration, and $u$ is jerk.
+ **Non-zero boundary states**. For example, we can require that the trajectory moves from $(x_1,x_2,x_3)=(1,-0.375,4)$ to $(x_1,x_2,x_3)=(-0.1,0.1,-1)$ in a 3rd-order problem.
+ **High computational efficiency**. The computation time of MIM can be significantly reduced compared to existing methods, where **3rd-order** problems require only about **0.2~0.8 ms** (in release mode), and the computation time for **4th-order** problems is reduced by **two orders of magnitude** compared to the existing optimization-based methods.
+ **High trajectory quality.** No chattering exists.
+ **High success rate**. 100% success rate for problems of order $n\leq4$.
+ **Strict/near time-optimality**. 100% optimality for problems of order $n\leq3$. 99.88% optimality for 4h-order problems.

## Statement

+ To achieve ultimate performance, I have rewritten MIM in **C++**. The paper used MATLAB, so this repository will even significantly outperform the effects claimed in the paper.
+ The development based on C++ is a massive undertaking, so I will proceed with a phased approach and gradually open-source these features. Noting that third-order trajectories are currently the most widely used, I have initially open-sourced the third-order trajectory. Higher-order trajectories, including fourth-order and above, will be open-sourced in the future as time permits. If you urgently need them, please contact me.
+ The project depends on: [GNU Scientific Library](https://www.gnu.org/software/gsl/).
+ If you use the MIM library for academic purposes, please cite our paper in IEEE TAC. If you have commercial needs, please contact me through wang-yn22@mails.tsinghua.edu.cn.

## Step-by-Step Guidance

### Step 1. Set the constraint and the boundary states

If you consider a jerk-limited (3rd-order) trajectory, please use:

```c++
#include "Planner.h"
// ...
int order = 3; // 3rd-order problem
double M_max[] = { 1.0, 1.0, 1.5, 4.0 }; // Maximal jerk, acceleration, velocity, and position
double M_min[] = { -1.0, -1.0, -1.5, -4.0 }; // Minimal jerk, acceleration, velocity, and position
Constraint constraint;
constraint.copy(order, M_max, M_min);
double x0[3] = { 1, -3.0 / 8.0, 4 }; // The initial acceleration, velocity, and position
double xf[3] = { -0.1, 0.1, -1 }; // The terminal acceleration, velocity, and position
```

If you need an infty, then code like this:

```c++
#include <limits>
// ...
double M_max[] = { 1.0, 1.0, 1.5, numeric_limits<double>::infinity() };
double M_min[] = { -1.0, -1.0, -1.5, -numeric_limits<double>::infinity() };
```

### Step 2.  Plan the trajectory

```c++
vector<arc> arcs = Planner::plan(order, x0, xf, constraint, true);
```

In `vector<arc> Planner::plan(int order, double* x0, double* xf, Constraint& constraint, bool flag_consider_position)`, 

+ `(int) order`: the order of the problem.
+ `(double*) x0`: the initial state vector with a length of `order`.
+ `(double*) xf`: the terminal state vector with a length of `order`.
+ `(Constraint&) constraint`: the constraint.
+ `(bool) flag_consider_position`: If it is `false`, then we ignore the position constraint, i.e., `M_max[order]` and `M_min[order]`  is considered as `infty`.

`arc` is a struct.

+ `arc.order` means that the constraint of which order is active. In 3rd-order problems, `arc.order==0` means that it is an unconstrained arc where the jerk is maximal or minimal. `arc.order==1` means that it is a constrained arc where the acceleration is maximal or minimal. `arc.order==2` means that it is a constrained arc where the velocity is maximal or minimal. 
+ `arc.sign==true`  if the control/state is maximal. `arc.sign==false`  if the control/state is minimal.
+ `arc.time` is the motion time of this arc.
+ `arc.tangent==0` is it is an arc.
+ `arc.tangent>0` is it is a tangent marker (of tangent order `arc.tangent`) instead of an arc. In this case, `arc.time==0` holds. The trajectory reaches the boundary here.

### Step 3. Interpolate the trajectory

```c++
double Ts = 1e-3; // the sample time
double T0 = 0; // 0<=T0<Ts is allowed.
Interpolator interpolator(order); // set the order
interpolator.interpolate(x0, arcs.data(), arcs.size(), constraint, Ts); // interpolate the trajectory
```

In this example of order 3, you can get:

+ `(vector<double>) interpolator.buffer[0]`: the control (jerk) at `{T0,T0+Ts,T0+2*Ts,...}`.
+ `(vector<double>) interpolator.buffer[1]`: the acceleration at `{T0,T0+Ts,T0+2*Ts,...}`.
+ `(vector<double>) interpolator.buffer[2]`: the velocity at `{T0,T0+Ts,T0+2*Ts,...}`.
+ `(vector<double>) interpolator.buffer[3]`: the position at `{T0,T0+Ts,T0+2*Ts,...}`.

If you want to output the trajectory as a csv file, you can type like this:

```c++
interpolator.write_csv(R"(..\data\3rd_order\0102010.csv)");
```

### Example 1. 3rd-order problems with 7 arcs.

```c++
#include "Planner.h"
#include "Interpolator.h"
using namespace std;

int main() {
	int order = 3;
	double M_max[4] = { 1.0, 1.0, 1.5, 4.0 };
	double M_min[4] = { -1.0, -1.0, -1.5, -4.0 };
	Constraint constraint;
	constraint.copy(order, M_max, M_min);
	double x0[3] = { 1, -3.0 / 8.0, 4 };
	double xf[3] = { -0.1, 0.1, -1 };

	vector<arc> arcs = Planner::plan(order, x0, xf, constraint, true);

	double Ts = 1e-3;
	Interpolator interpolator(order);
	interpolator.interpolate(x0, arcs.data(), arcs.size(), constraint, Ts);
	interpolator.write_csv(R"(..\data\3rd_order\0102010.csv)");
}
```

Then, you can get the following results: (`arc` is a struct)

| arc       | order | sign    | tangent | time    |
| --------- | ----- | ------- | ------- | ------- |
| `arcs[0]` | 0     | `false` | 0       | 2 |
| `arcs[1]` | 1     | `false`  | 0     | 0.625    |
| `arcs[2]` | 0     | `true` | 0       | 1 |
| `arcs[3]` | 2     | `false`  | 0       | 0.897994 |
| `arcs[4]` | 0     | `true` | 0       | 1 |
| `arcs[5]` | 1     | `true`  | 0       | 0.605 |
| `arcs[6]` | 0     | `false`  | 0       | 1.1 |

That means the augmented switching law is

$$
\underline{01}\overline{0}\underline{2}\overline{01}\underline{0}
$$

+ Jerk is `M_min[0]` during $0<t<2$.
+ Acceleration is `M_min[1]` and jerk is `0` during $2<t<2.625$.
+ Jerk is `M_max[0]` during $2.625<t<3.625$.
+ Velocity is `M_min[2]`, acceleration is `0`, and jerk is `0` during $3.625<t<4.522994$.
+ Jerk is `M_max[0]` during $4.522994<t<5.522994$.
+ Acceleration is `M_max[1]` and jerk is `0` during $5.522994<t<6.127994$.
+ Jerk is `M_min[0]` during $6.127994<t<7.227994$.

<p align="center">
	<img src="https://github.com/user-attachments/assets/8016f5ef-3003-4b1f-9ec3-e8860687e358" alt="Example_1" width="500" />
</p>
<p align="center">
	<b>Figure.</b> The trajectory in Example 1.
</p>

### Example 2. 3rd-order problems with 5 arcs and 1 tangent marker.

```c++
#include "Planner.h"
#include "Interpolator.h"
using namespace std;

int main() {
	int order = 3;
	double M_max[4] = { 1.0, 1.0, 1.5, 4.0 };
	double M_min[4] = { -1.0, -1.0, -1.5, -4.0 };
	Constraint constraint;
	constraint.copy(order, M_max, M_min);
	double x0[3] = { 1,-0.375, 3.999 };
	double xf[3] = { 0, 0, 4 };

	vector<arc> arcs = Planner::plan(order, x0, xf, constraint, true);

	double Ts = 1e-3;
	Interpolator interpolator(order);
	interpolator.interpolate(x0, arcs.data(), arcs.size(), constraint, Ts);
	interpolator.write_csv(R"(..\data\3rd_order\00_3_2_000.csv)");
}
```

Then, you can get the following results: (`arc` is a struct)

| arc       | order | sign    | tangent | time    |
| --------- | ----- | ------- | ------- | ------- |
| `arcs[0]` | 0     | `false` | 0       | 1.38027 |
| `arcs[1]` | 0     | `true`  | 0       | 0.18226 |
| `arcs[2]` | 3     | `false` | 2       |         |
| `arcs[3]` | 0     | `true`  | 0       | 0.39503 |
| `arcs[4]` | 0     | `false` | 0       | 0.33565 |
| `arcs[5]` | 0     | `true`  | 0       | 0.13862 |

That means the augmented switching law is

$$
\underline{0}\overline{0}(\overline{3},2)\overline{0}\underline{0}\overline{0}
$$

+ Jerk is `M_min[0]` during $0<t<1.38027$.
+ Jerk is `M_max[0]` during $1.38027<t<1.56253$.
+ The position is tangent to `M_max[4]` of order `2`. At $t=1.56253$, the position is `M_max[4]`, the velocity is `0`, and the acceleration is `<0`.
+ Jerk is `M_max[0]` during $1.56253<t<1.95756$.
+ Jerk is `M_min[0]` during $1.95756<t<2.29321$.
+ Jerk is `M_max[0]` during $2.29321<t<2.43183$.

<p align="center">
	<img src="https://github.com/user-attachments/assets/2d919ba0-1895-46dc-9aab-fd0f5f534a90" alt="Example_2" width="500" />
</p>
<p align="center">
	<b>Figure.</b> The trajectory in Example 2.
</p>

