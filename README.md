# Quaternion Safety Filter

MATLAB implementation of the robust safety filter for nonlinear quaternion attitude dynamics presented in:

> **Robust Safety Filter Synthesis for Quaternion Attitude Dynamics via LMI-Based Ellipsoidal Invariant Sets**  
> Reza Pordal, Alireza Sharifi, Ali Baniasad  
> Department of Aerospace Engineering, Sharif University of Technology

## Overview

This repository provides the code to reproduce all simulation results in the paper. The framework synthesizes a maximal ellipsoidal robust controlled invariant (RCI) set and an associated state-feedback backup controller for quaternion attitude dynamics by solving a single convex semidefinite program (SDP). Exact closed-form sector bounds on the quaternion kinematic nonlinearity are derived analytically and embedded into the LMI via the S-procedure, providing formal safety guarantees for the full nonlinear system without linearization.

## Requirements

- MATLAB R2020b or later
- [CVX](http://cvxr.com/cvx/) (convex optimization toolbox)
- [MOSEK](https://www.mosek.com/) solver (recommended; free academic license available)

## File Structure

| File | Description |
|------|-------------|
| `quadrotor_example_01.m` | **Scenario I**: set-point tracking, small goal at (2, 1) m — nominal operation |
| `quadrotor_example_02.m` | **Scenario II**: set-point tracking, distant goal at (25, 15) m — large initial errors triggering filter |
| `quadrotor_example_03.m` | **Scenario III**: circular trajectory (radius 5 m, 1 rad/s) — persistent attitude excitation |
| `loadSystemParameters.m` | Crazyflie 2.0 quadrotor parameters (mass, inertia) |
| `plotSafetyResults.m` | 2×3 paper figure: position, Euler angles, torques, quaternion norm, Lyapunov function, sector bound ratio |
| `formatFigureIEEE.m` | Figure formatting utility for IEEE-style plots |

## Running the Simulations

Open MATLAB, navigate to this folder, and run any of the scenario scripts:

```matlab
% Scenario I — small setpoint
quadrotor_example_01

% Scenario II — large setpoint (demonstrates filter intervention)
quadrotor_example_02

% Scenario III — circular trajectory
quadrotor_example_03
```

Each script:
1. Constructs the linearized quaternion attitude dynamics for the Crazyflie 2.0
2. Computes the exact sector bound `gamma` from the admissible rotation angle
3. Solves the LMI optimization (via CVX + MOSEK) to find the RCI set and backup controller
4. Runs the nonlinear closed-loop simulation under bounded disturbances
5. Plots the results

## System Parameters

The quadrotor model is a Crazyflie 2.0 nano quadrotor:

| Parameter | Value |
|-----------|-------|
| Mass | 0.028 kg |
| Inertia (Ixx, Iyy) | 16.6 × 10⁻⁶ kg·m² |
| Inertia (Izz) | 29.3 × 10⁻⁶ kg·m² |
| Max rotation angle | 40° |
| Max angular rate | 5 rad/s |
| Max torque per axis | 1 × 10⁻⁴ N·m |
| Disturbance bound | 1 × 10⁻⁵ N·m |

## Citation

If you use this code, please cite:

```bibtex
@inproceedings{pordal2026quaternion,
  title     = {Robust Safety Filter Synthesis for Quaternion Attitude Dynamics via {LMI}-Based Ellipsoidal Invariant Sets},
  author    = {Pordal, Reza and Sharifi, Alireza and Baniasad, Ali},
  booktitle = {Proceedings of the IEEE Conference},
  year      = {2026}
}
```

## License

MIT License — see [LICENSE](LICENSE) for details.
