# Simulated Annealing 

Simulated Annealing (SA) is a stochastic optimization technique inspired by the physical annealing process in solids. It is especially effective for black-box optimization problems where the objective function is costly to evaluate, non-differentiable, or unknown in closed form.

The algorithm simulates how a material cools and gradually settles into a low-energy state. In this analogy:

- The **state space** corresponds to all possible solutions.
- The **objective function** represents the system’s energy, which the algorithm seeks to minimize.
- The **temperature** parameter controls the probability of accepting worse solutions to avoid getting trapped in local minima.

At each iteration, a new candidate solution is generated. If it improves the objective function, it is accepted unconditionally. Otherwise, it is accepted with probability:

$$ \
P(\text{accept}) = \exp\left(-\frac{\Delta f}{k \cdot T}\right)
\ $$

where:

- $$\(\Delta f\)$$ is the increase in the objective function value,
- $$\(k\)$$ is a constant analogous to Boltzmann’s constant,
- $$\(T\)$$ is the current temperature.

This acceptance rule, known as the *Metropolis criterion*, allows the algorithm to balance exploration and exploitation, gradually focusing the search as the temperature decreases.


## Multithreading version (OMP) 
### Main
> ⚙️ **Executable Overview** – This section explains how to run the main OpenMP Simulated Annealing implementation.


#### Execution
An example invocation of the `main` executable for the OpenMP Simulated Annealing implementation (all the detailed flag for main execution are in the main readme page):

```bash
./build/main -a simulated_omp -d 2 -i 100 -f sphere
```

## 🔍 Convergence Tests

This section documents a set of tests used to verify the **convergence behavior of the Simulated Annealing algorithm** on various classical objective functions. The tests are implemented using GoogleTest and check both the **accuracy of the solution** and the **closeness of the estimated minimum position** to the known global optimum.

Each test asserts that the algorithm achieves a result **sufficiently close to the known global minimum**, with different error thresholds tailored to the specific characteristics of each objective function.

### ⚠️ Parameter Sensitivity

The **algorithm parameters** (such as initial temperature, cooling rate, step size, etc.) are **not one-size-fits-all**—they are **tuned per function**, because each objective function presents unique challenges, such as:

- **Local curvature**: For instance, the Rastrigin function has many sharp local minima, while the Sphere function is smooth and convex.
- **Problem conditioning**: The Rosenbrock function has a narrow, curved valley that is hard to explore without fine-grained control.
- **Search domain scale** and **number of local minima**: Parameters like `step_size` and `temperature_scale` greatly affect the algorithm's ability to escape local optima and efficiently explore the space.

Choosing good parameters requires balancing exploration (via temperature and step size) with convergence (via cooling and step decay), and often involves empirical tuning.

---

#### ✅ Functions Tested

| Test Name                      | Objective Function    | Convergence Criterion                          |
|-------------------------------|------------------------|------------------------------------------------|
| `SaConvergence.Sphere`        | Sphere                 | ‖f(x) – 0‖ ≤ 1 × 10⁻³ <br>‖x – 0‖ ≤ 5 × 10⁻²     |
| `SaConvergence.EuclideanDistance` | Euclidean Distance  | ‖f(x) – 0‖ ≤ 3 × 10⁻³ <br>‖x – 0‖ ≤ 3 × 10⁻³     |
| `SaConvergence.Rosenbrock`    | Rosenbrock             | ‖f(x) – 0‖ ≤ 1 × 10⁻³ <br>‖x – x*‖ ≤ 1 × 10⁻²   |
| `SaConvergence.Rastrigin`     | Rastrigin              | ‖f(x) – 0‖ ≤ 1 × 10⁻² <br>‖x – 0‖ ≤ 2 × 10⁻²     |


---

### 🛠️ Test Setup

Each test shares the same general structure but uses function-specific parameters. All tests invoke the core routine:

```cpp
algorithm::run_simulated_annealing(
    dimensions,
    max_iterations,
    dwell,
    initial_temperature,
    temperature_scale,
    initial_step_size,
    step_size_scale,
    boltzmann_k,
    initial_guess,
    lower_bound,
    upper_bound,
    objective_function,
    seed,
    n_threads,
    true  // verbose
);

#### Test result


#### Summary of test duration


### Benchmark
