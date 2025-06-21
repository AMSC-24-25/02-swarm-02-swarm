# Firefly + BFGS Optimization with CPU & GPU Acceleration

This project implements a hybrid metaheuristic optimization algorithm that combines the Firefly Algorithm (FA) with the BFGS (quasi-Newton method) for high-precision refinement

#### Firefly Algorithm Update Rule

$$
x_i(t+1) = x_i(t) + \beta_0 e^{-\gamma r_{ij}^2} (x_j(t) - x_i(t)) + \alpha \cdot \varepsilon
$$
where: 
$$
\begin{aligned}
x_i(t) &\quad \text{Position of firefly } i \text{ at iteration } t \\
x_j(t) &\quad \text{Position of a brighter firefly } j \\
\beta_0 &\quad \text{Attractiveness at distance } 0 \\
\gamma &\quad \text{Light absorption coefficient} \\
r_{ij} &\quad \text{Euclidean distance between } x_i \text{ and } x_j \\
\alpha &\quad \text{Randomness scaling factor} \\
\varepsilon &\quad \text{Random vector from } \mathcal{U}(-0.5, 0.5)
\end{aligned}
$$

BFGS refines the best solution found by FA using a gradient-based update:

$$
x_{k+1} = x_k - H_k \nabla f(x_k)
$$

Where:

$$
\begin{aligned}
x_k &\quad \text{Current solution at iteration } k \\
\nabla f(x_k) &\quad \text{Gradient of the objective function at } x_k \\
H_k &\quad \text{Approximate inverse Hessian matrix} \\
x_{k+1} &\quad \text{Updated (refined) solution}
\end{aligned}
$$

---


🔧 CMake Build

To build with CUDA support (by dafault Cuda is OFF):
```bash
cmake -S . -B build -DENABLE_CUDA=ON
```

### Example Execution
```bash
./build/main -a firefly_bfgs -d 2 -n 100 -i 500 -f sphere
```
---
📄 **Full explanation and benchmarks available in** [`reportFireFlyBFGS.pdf`](https://github.com/AMSC-24-25/02-swarm-02-swarm/blob/main/doc/reportFireFlyBFGS.pdf)


