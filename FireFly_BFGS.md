# Firefly + BFGS Optimization with CPU & GPU Acceleration

This project implements a hybrid metaheuristic optimization algorithm that combines the Firefly Algorithm (FA) with the BFGS (quasi-Newton method) for high-precision refinement

x_i(t+1) = x_i(t) + β₀ * exp(-γ * r_ij²) * (x_j(t) - x_i(t)) + α * rand()
Where:

β₀ = base attractiveness,

γ = light absorption coefficient,

r_ij = distance between fireflies i and j,

α = randomness factor,

rand() = random vector ∈ [-0.5, 0.5].

BFGS refines the best solution found by FA using a gradient-based update:
x_{k+1} = x_k - H_k * ∇f(x_k)

Where:

H_k = approximation of the inverse Hessian,

∇f(x_k) = gradient at point x_k.


🔧 CMake Build
To build with CUDA support (dafault is OFF):
```bash
cmake .. -DENABLE_CUDA=ON
```

### Example Execution
```bash
./build/main -a firefly_bfgs -d 2 -n 100 -i 500 -f sphere
```

📄 Documentation
📘 For a full explanation of the algorithm, design choices, and benchmarks, see doc/reportFireFlyBFGS.pdf.

