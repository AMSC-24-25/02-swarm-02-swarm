# Firefly + BFGS Optimization with CPU & GPU Acceleration

This project implements a hybrid metaheuristic optimization algorithm that combines the Firefly Algorithm (FA) with the BFGS (quasi-Newton method) for high-precision refinement

#### Firefly Algorithm Update Rule
```latex
x_i(t+1) = x_i(t) + \beta_0 e^{-\gamma r_{ij}^2} (x_j(t) - x_i(t)) + \alpha \cdot \varepsilon
```

Where:
- `β₀`: base attractiveness  
- `γ`: light absorption coefficient  
- `r_ij`: Euclidean distance between fireflies `i` and `j`  
- `α`: randomness factor  
- `rand()`: random vector in [-0.5, 0.5]

BFGS refines the best solution found by FA using a gradient-based update:

```latex
x_{k+1} = x_k - H_k \nabla f(x_k)

```

Where:
- `H_k`: approximation of the inverse Hessian  
- `∇f(x_k)`: gradient at point `x_k`


---


🔧 CMake Build
To build with CUDA support (dafault is OFF):
```bash
cmake .. -DENABLE_CUDA=ON
```

### Example Execution
```bash
./build/main -a firefly_bfgs -d 2 -n 100 -i 500 -f sphere
```
---
📄 **Full explanation and benchmarks available in** [`reportFireFlyBFGS.pdf`](https://github.com/AMSC-24-25/02-swarm-02-swarm/blob/main/doc/reportFireFlyBFGS.pdf)


