# Manifold-Intercept Method

[![License](https://img.shields.io/badge/License-Boost_1.0-lightblue.svg)](https://www.boost.org/LICENSE_1_0.txt)

**Coming soon...**

The Manifold-Intercept Method (MIM) is a one-axis point-to-point trajectory planning method with full box state constraints. The inistial state and the terminal state are allow to be non-zero.

MIM can achieve time-optimalilty for 3rd order problems, i.e., jerk-limited trajectories. It can plan suboptimal trajectory under full box state constraints for 4th order and higher-order problems. For 4th order problems, i.e., snap-limited trajectories, the relative error on terminal time between MIM and the optimal solution with chattering is usually less than 0.2% in practice; if apply MIM to substitude the chattering period in time-optimal 4th order trajectories, then the relative error on terminal time between MIM and the optimal solution will be less than 0.12% strictly. 

The MIM method is based on the following paper. It is under reviewed now. Therefore, it will be open-source once the paper is public.

```
@article{wang2024time,
  title={Time-Optimal Control for High-Order Chain-of-Integrators Systems with Full State Constraints and Arbitrary Terminal States (Part I, Extended Version)},
  author={Wang, Yunan and Hu, Chuxiong and Li, Zeyang and Lin, Shize and He, Suqin and Wang, Ze and Zhu, Yu},
  journal={arXiv preprint arXiv:2311.07039},
  year={2024}
}
```
