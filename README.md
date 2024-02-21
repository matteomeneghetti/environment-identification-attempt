# Environment Identification attempt

This repository summarizes an attempt in system identification that I made during my time at the Altair Robotics Laboratory.

Specifically, it is about the synthesis of an estimator able to evaluate both the tightly coupled environment stiffness and equilibrium position. This has proven to be quite challenging through the use of linear estimators.

The environment is modeled through the use of a MATLAB S-function as a second order system $\tau_e = K_e(\theta_e - \theta_0) + b_e\dot{\theta} + J_e\ddot{\theta}$ where $K_e$, $b_e$, $J_e$ are respectively the environment stiffness, damping and moment of inertia and $\theta_0$ is the equilibrium position of the ideal spring.

Since we are mostly concerned about an environment which continously changes $K_e$ and $\theta_0$, both $b_e$ and $J_e$ can be kept equal to zero, leaving only

$$ \tau_e = K_e(\theta_e - \theta_0) $$

A basic implementation of a linear Recursive Least-Square regressor with forgetting factor $\lambda$ is provided, which has been empirically proven equal to the one offered by the MATLAB System Identification Toolbox.
