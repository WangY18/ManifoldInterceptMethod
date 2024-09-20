# Manifold-Intercept Method

[![License](https://img.shields.io/github/license/WangY18/ManifoldInterceptMethod)](https://www.gnu.org/licenses/gpl-3.0.en.html)

**Coming soon...**

The Manifold-Intercept Method (MIM) is a one-axis point-to-point trajectory planning method with full box state constraints. The inistial state and the terminal state are allow to be non-zero.

MIM can achieve time-optimalilty for 3rd order problems, i.e., jerk-limited trajectories. It can plan suboptimal trajectory under full box state constraints for 4th order and higher-order problems. For 4th order problems, i.e., snap-limited trajectories, the relative error on terminal time between MIM and the optimal solution with chattering is usually less than 0.2% in practice; if apply MIM to substitude the chattering period in time-optimal 4th order trajectories, then the relative error on terminal time between MIM and the optimal solution will be less than 0.12% strictly. 

The MIM method is based on the following paper. I am working on rewrite the MIM to make it open-source now. I believe that it can be open-source in 2024.11.

```
@article{wang2024time,
  title={Time-Optimal Control for High-Order Chain-of-Integrators Systems with Full State Constraints and Arbitrary Terminal States},
  author={Wang, Yunan and Hu, Chuxiong and Li, Zeyang and Lin, Shize and He, Suqin and Zhu, Yu},
  journal={IEEE Transactions on Automatic Control},
  year={2024}
}
```
