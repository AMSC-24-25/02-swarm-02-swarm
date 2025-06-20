# Simulated Annealing 

Simulated Annealing (SA) is a stochastic optimization technique inspired by the physical annealing process in solids. It is especially effective for black-box optimization problems where the objective function is costly to evaluate, non-differentiable, or unknown in closed form.

The algorithm simulates how a material cools and gradually settles into a low-energy state. In this analogy:

- The **state space** corresponds to all possible solutions.
- The **objective function** represents the system’s energy, which the algorithm seeks to minimize.
- The **temperature** parameter controls the probability of accepting worse solutions to avoid getting trapped in local minima.

At each iteration, a new candidate solution is generated. If it improves the objective function, it is accepted unconditionally. Otherwise, it is accepted with probability:

$$ \[
P(\text{accept}) = \exp\left(-\frac{\Delta f}{k \cdot T}\right)
\] $$

where:

- \(\Delta f\) is the increase in the objective function value,
- \(k\) is a constant analogous to Boltzmann’s constant,
- \(T\) is the current temperature.

This acceptance rule, known as the *Metropolis criterion*, allows the algorithm to balance exploration and exploitation, gradually focusing the search as the temperature decreases.


## Multithreading version (OMP) 
### Main
> ⚙️ **Executable Overview** – This section explains how to run the main OpenMP Simulated Annealing implementation.


#### Execution
An example invocation of the `main` executable for the OpenMP Simulated Annealing implementation (all the detailed flag for main execution are in the main readme page):

```bash
./build/main -a simulated_omp -d 2 -i 100 -f sphere
```
