# Simulated Annealing 
Simulated Annealing (SA) is a stochastic optimization method inspired by the physical process of annealing in solids. It is particularly useful for black-box optimization problems where the objective function is expensive, non-differentiable, or unknown in closed form.

The algorithm mimics how a material cools down and settles into a low-energy state. In this context:

The state space represents all possible solutions.

The objective function plays the role of the system's energy, which we aim to minimize.

The temperature controls the likelihood of accepting worse solutions to escape local minima.

At each step, a new solution is proposed. If it improves the objective, it is accepted. Otherwise, it's accepted with probability:

P(accept) = exp( -Δf / (k * T) )
Where:

Δf is the increase in objective value,

k is a constant (analogous to Boltzmann’s constant),

T is the current temperature.

This strategy, known as the Metropolis criterion, helps balance exploration and exploitation, gradually focusing the search as the temperature decreases.

## Multithreading version (OMP) 
### Main
> ⚙️ **Executable Overview** – This section explains how to run the main OpenMP Simulated Annealing implementation.


#### Execution
An example invocation of the `main` executable for the OpenMP Simulated Annealing implementation (all the detailed flag for main execution are in the main readme page):

```bash
./build/main -a differential_omp -d 2 -n 100 -i 100 -f sphere
```