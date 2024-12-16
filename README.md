# SwarmSearch
Particle swarm optimization (PSO) is a stochastic population-based derivative-free optimization method which offers a high parallelization.

The minimization problem can be defined as:

$$
min f(x), f : \R^n \rightarrow \R
$$

such that

$$
x \in \R^n, g_i(x) ≤ 0, i = {1...M}
$$

## The algorithm
The algorithm involves creating a swarm of particles which, in turn, represent various solutions moving in the n-dimensional input space of the objective function. At each iteration, all the particles' velocities and positions are updated with the following formulas:

$$
v_i^{t+1} = \omega v_i^t + c_1 r_1 \left( p_{\text{best},i}^t - x_i^t \right) + c_2 r_2 \left( g_{\text{best}}^t - x_i^t \right)
$$

$$
x_i^{t+1} = x_i^t + v_i^{t+1}
$$

The global minimum coincides with the particle with the "best" value across the swarm.

## The implementation
The main classes involved in this project are:
- Swarm: keeps track of the best positions of all the particles and updates them accordingly. It also keeps the main coefficients of the algorithm.
- Particle: just contains its position and velocity, along with the best position across iterations, and is updated by the swarm in each iteration.
- ObjectiveFunction: abstract class to represent a function to be minimized.

The main area in which we parallelized was the `Swarm::updateParticles()` method which is embarrassingly parallel, meaning that a single OpenMP `#pragma` directive was sufficient to achieve near-linear speedup. Another area involved in the optimization was the `Swarm::findBestFitness()` method which performs a reduction across all particles in the swarm to find the best position to be used in the successive iterations. This kind of reduction was not so straightforward since it needed to find only the index of the best particle (instead of its value), so we implemented it "manually" through a parallel for followed by an update of the global result inside a `critical` section.

To check the algorithm's correctness and performance, we tested it against different functions with increasingly harder to minimize:
 - Sphere: the easiest one, a single global minimum.
 - EuclideanDistance: a single global minimum.
 - [Rosenbrock](https://en.wikipedia.org/wiki/Rosenbrock_function): the "banana" function, notoriously hard to find the minimum inside the main valley.
 - [Rastrigin](https://en.wikipedia.org/wiki/Rastrigin_function): a lot of local minima clustered around the central global minimum.

## How to compile and run
Compile with:
```bash
make
```

Run with:
```bash
./build/main
```

Check that everything works correctly
```bash
make test
```
