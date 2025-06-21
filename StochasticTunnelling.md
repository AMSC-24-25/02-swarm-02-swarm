#  Stochastic Tunnelling Approach
Stochastic Tunnelling (ST) is an optimization algorithm designed for determining the global minima of complex and rugged energy landscapes. It is a generic physically motivated generalization of
simulated annealing. This approach circumvents the freezing problem which arises when the energy difference between “adjacent” local minima on the energy surface is much smaller than the energy of 
intervening transition states separating them.The physical idea behind the stochastic tunneling method is to allow the particle to “tunnel” forbidden regions, once it has been determined
that they are irrelevant for the low-energy properties of the problem. This can be accomplished by applying the transformation:

$$
f_{\text{STUN}}(x) = 1 - \exp[-\gamma(f(x) - f_0)],
$$

where $f_0$ is the lowest minimum encountered thus far. The effective potential preserves the locations of all minima, but maps the entire energy space from $f_0$ to the maximum
of the potential onto the interval [0, 1]. The degree of steepness of the cutoff of the high-energy regions is controlled by the tunneling parameter $\gamma$. To illustrate the physical content of the transformation we consider a Monte Carlo process at some fixed inverse
temperature $\beta$. A step from $x_1$ to $x_2$ with $\Delta = f(x_1) - f(x_2)$ ­  is accepted with probability $\tilde{w}_{1\to 2}=\exp(-\tilde{\beta}\Delta)$. In the process the temperature rises rapidly when the local energy is larger than $f_0$ and the particle diffuses (or tunnels) freely through potential barriers of arbitrary height. As better and better minima are found, ever larger
portions of the high-energy part of the function surface are flattened out. In analogy to the simulated annealing approach this behavior can be interpreted as a self-adjusting cooling schedule that is optimized as the simulation proceeds. Parameter $\beta$ is adjusted during the simulation. If a short-time moving average of the $f_{stun}$ exceeds the threshold $f_{\text{thresh}}$, $\beta$ is reduced by some fixed factor, otherwise it is increased. \beta is updated in the following manner:




## Multithreading version (OMP)
### Main
> ⚙️ **Executable Overview** – This section explains how to run the main OpenMP Differential Evolution implementation.


#### Execution
An example invocation of the `main` executable for the OpenMP Differential Evolution implementation (all the detailed flag for main execution are in the main readme page):

```bash
./build/main -a differential_omp -d 2 -n 100 -i 100 -f sphere
```

#### Sample output
```console
Iteration n. 1 / 100
  Current minimum:
  f(-7.798279e-01, 2.117076e+00) = 5.090143e+00

Iteration n. 2 / 100
  Current minimum:
  f(-7.798279e-01, 2.117076e+00) = 5.090143e+00

…

Iteration n. 99 / 100
  Current minimum:
  f(1.976930e-13, 1.210598e-13) = 5.373799e-26

Iteration n. 100 / 100
  Current minimum:
  f(1.976930e-13, 1.210598e-13) = 5.373799e-26

Minimum found:
  f(1.976930e-13, 1.210598e-13) = 5.373799e-26
Total execution time: 0.020753 seconds
```
#### Comments

- **Progress reporting**: Each iteration prints the current best solution (coordinates and objective value).  
- **Convergence behavior**: The objective value decreases rapidly, reaching ~5×10⁻²⁶ by iteration 100.  
- **Performance**: 100 iterations on a 2-D Sphere function with 100 candidates and 5 threads complete in ~0.02 s.

### Test
> ✅ **Test Suite Summary** – Automated tests  to validate correctness and convergence.


A comprehensive test suite is provided to verify the correctness and convergence properties of the Differential Evolution implementation. All tests are written in C++17 using the GoogleTest framework and exercise the algorithm on four classic benchmark functions.

#### Function tested 
| Test Name                         | Objective Function | Convergence Criterion                |
| --------------------------------- | ------------------ | ------------------------------------ |
| `DeConvergence.Sphere`            | Sphere             | ‖f(x) – 0‖ ≤ 1 × 10⁻³                |
| `DeConvergence.EuclideanDistance` | Euclidean Distance | ‖f(x) – 0‖ ≤ 1 × 10⁻³                |
| `DeConvergence.Rosenbrock`        | Rosenbrock         | Position error ≤ 0.1 (narrow valley) |
| `DeConvergence.Rastrigin`         | Rastrigin          | ‖f(x) – 0‖ ≤ 1 × 10⁻³                |


#### Test Setup

- **Framework:** GoogleTest  
- **Dimensions:** 2  
- **Population size:** 100 candidates  
- **Max iterations:** 1 000  
- **Random seed:** 42  
- **Search bounds:** \[-10, 10\]  
- **Differential weight (F):** 0.5  
- **Crossover rate (CR):** 0.8  
- **Threads:** 5  

Each test invokes:
```cpp
auto result = algorithm::run_differential_evolution(
    dimensions,
    num_candidates,
    lower_bound,
    upper_bound,
    seed,
    max_iterations,
    F,
    CR,
    objectiveFunction,
    /* num_threads= */ 5,
    /* verbose= */ false
);
```
#### Test result 
All four tests passed:

```console
[==========] Running 4 tests from 1 test suite.
[----------] Global test environment set-up.
[ RUN      ] DeConvergence.Sphere
[       OK ] DeConvergence.Sphere (81 ms)
[ RUN      ] DeConvergence.EuclideanDistance
[       OK ] DeConvergence.EuclideanDistance (67 ms)
[ RUN      ] DeConvergence.Rosenbrock
[       OK ] DeConvergence.Rosenbrock (65 ms)
[ RUN      ] DeConvergence.Rastrigin
[       OK ] DeConvergence.Rastrigin (66 ms)
[----------] 4 tests from DeConvergence (280 ms total)
[  PASSED  ] 4 tests.
```


#### Summary of test duration
| Test                   | Duration |
| ---------------------- | -------: |
| Sphere                 |    81 ms |
| EuclideanDistance      |    67 ms |
| Rosenbrock             |    65 ms |
| Rastrigin              |    66 ms |
| **Total elapsed time** |   280 ms |


### Benchmark
> 📊 **Performance Benchmarking** – Evaluate speedup and scaling on different workloads and thread counts.

The algorithm was tested with varying numbers of threads and different population sizes (called "creatures").

####  Time vs Number of Creatures

![](benchResult/figures/time_vs_creatures.png)

This plot shows how the execution time (in seconds, log scale) varies as the number of creatures increases, for 1, 8, and 16 threads:

- With **1 thread**, the time grows consistently with the number of creatures.
- With **8 and 16 threads**, the performance improves significantly when the number of creatures is sufficiently large.
- For **small numbers of creatures**, multithreading introduces overhead that can outweigh the benefits of parallelization.
- The curves flatten in the low-creature region for 8 and 16 threads, showing that overhead dominates until the workload becomes heavy enough.

####  Strong Speedup vs Number of Threads

![](benchResult/figures/speedup_vs_threads.png)

This plot shows **strong scaling**, i.e., how speedup changes as the number of threads increases, for a fixed number of creatures (4, 512, and 1024):

- **4 creatures:** Speedup decreases with more threads due to overhead dominating such a small workload.
- **512 creatures:** Modest speedup is observed, showing that parallelism becomes more effective.
- **1024 creatures:** Best performance scaling, achieving a speedup > 2.5× with 16 threads.

####  Conclusions

- Parallel execution is effective only when the number of creatures is large enough to amortize thread management overhead.
- For small problem sizes, single-threaded execution remains more efficient.
- The algorithm demonstrates **scalability potential**, especially on larger workloads.


## Multiprocessing version (MPI)

> **Note:**
> Follow the instructions in the main readme page to compile properly to use MPI.
> This MPI-based implementation is provided for educational purposes. In practice, MPI shines on distributed systems or clusters, whereas on a single multi-core machine it can incur additional communication overhead. As a result, we do not include a full benchmark suite here (and the Google Benchmark integration currently causes conflicts). What remains is:
> 1. an MPI-enabled `main` that runs Differential Evolution across multiple processes, and  
> 2. the same convergence tests you’ve already seen (using GoogleTest).

The algorithm is partitioned so that each MPI rank evolves its own subset of the population. At each generation:
1. every rank runs one DE update on its local candidates,  
2. each rank computes its best local solution,  
3. an `MPI_Allreduce` (with `MPI_MINLOC`) finds the global-best fitness and the rank that holds it,  
4. that rank broadcasts its best candidate vector to all others,  
5. the next generation proceeds using that shared global best.  

This design minimizes inter-process communication (only one `Allreduce` and one `Bcast` per iteration), at the cost of never mixing candidates between ranks.

### Main  
You can run it like so (full detail flag usage in main readme page):
```bash
mpirun -n 4 ./build/main -a differential_mpi -d 2 -n 100 -i 100 -f sphere -j 1
```
Sample output (only a few iterations shown):
```console
Iteration n. 1 / 100
  Current minimum:
  f(-8.384779e+00, 2.789538e+00) = 7.808604e+01

Iteration n. 2 / 100
  Current minimum:
  f(-7.100683e+00, 5.125734e+00) = 7.669284e+01

…

Iteration n. 98 / 100
  Current minimum:
  f(-1.217977e-12, -9.401532e-13) = 2.367357e-24

Iteration n. 99 / 100
  Current minimum:
  f(-1.217977e-12, -9.401532e-13) = 2.367357e-24

Iteration n. 100 / 100
  Current minimum:
  f(2.020037e-13, -1.754170e-13) = 7.157664e-26

Minimum found:
  f(2.020037e-13, -1.754170e-13) = 7.157664e-26  
Total execution time: 0.001724 seconds
```
### Test
You can run the tests like so:
```bash
mpirun -n 4 test_de_convergence_mpi
```

All the test are passed even in the multiprocessing version of the algorithm:
```console
[==========] Running 4 tests from 1 test suite.
[----------] Global test environment set-up.
[----------] 4 tests from DeConvergenceMPI
[ RUN      ] DeConvergenceMPI.Sphere
[==========] Running 4 tests from 1 test suite.
[----------] Global test environment set-up.
[----------] 4 tests from DeConvergenceMPI
[ RUN      ] DeConvergenceMPI.Sphere
[==========] Running 4 tests from 1 test suite.
[----------] Global test environment set-up.
[----------] 4 tests from DeConvergenceMPI
[ RUN      ] DeConvergenceMPI.Sphere
[==========] Running 4 tests from 1 test suite.
[----------] Global test environment set-up.
[----------] 4 tests from DeConvergenceMPI
[ RUN      ] DeConvergenceMPI.Sphere
[       OK ] DeConvergenceMPI.Sphere (8 ms)
[ RUN      ] DeConvergenceMPI.EuclideanDistance
[       OK ] DeConvergenceMPI.Sphere (8 ms)
[ RUN      ] DeConvergenceMPI.EuclideanDistance
[       OK ] DeConvergenceMPI.Sphere (8 ms)
[ RUN      ] DeConvergenceMPI.EuclideanDistance
[       OK ] DeConvergenceMPI.Sphere (8 ms)
[ RUN      ] DeConvergenceMPI.EuclideanDistance
[       OK ] DeConvergenceMPI.EuclideanDistance (9 ms)
[ RUN      ] DeConvergenceMPI.Rosenbrock
[       OK ] DeConvergenceMPI.EuclideanDistance (9 ms)
[ RUN      ] DeConvergenceMPI.Rosenbrock
[       OK ] DeConvergenceMPI.EuclideanDistance (9 ms)
[ RUN      ] DeConvergenceMPI.Rosenbrock
[       OK ] DeConvergenceMPI.EuclideanDistance (9 ms)
[ RUN      ] DeConvergenceMPI.Rosenbrock
[       OK ] DeConvergenceMPI.Rosenbrock (5 ms)
[ RUN      ] DeConvergenceMPI.Rastrigin
[       OK ] DeConvergenceMPI.Rosenbrock (5 ms)
[ RUN      ] DeConvergenceMPI.Rastrigin
[       OK ] DeConvergenceMPI.Rosenbrock (5 ms)
[ RUN      ] DeConvergenceMPI.Rastrigin
[       OK ] DeConvergenceMPI.Rosenbrock (5 ms)
[ RUN      ] DeConvergenceMPI.Rastrigin
[       OK ] DeConvergenceMPI.Rastrigin (5 ms)
[----------] 4 tests from DeConvergenceMPI (28 ms total)

[----------] Global test environment tear-down
[==========] 4 tests from 1 test suite ran. (28 ms total)
[  PASSED  ] 4 tests.
[       OK ] DeConvergenceMPI.Rastrigin (5 ms)
[----------] 4 tests from DeConvergenceMPI (28 ms total)

[----------] Global test environment tear-down
[==========] 4 tests from 1 test suite ran. (28 ms total)
[  PASSED  ] 4 tests.
[       OK ] DeConvergenceMPI.Rastrigin (5 ms)
[----------] 4 tests from DeConvergenceMPI (28 ms total)

[----------] Global test environment tear-down
[==========] 4 tests from 1 test suite ran. (28 ms total)
[  PASSED  ] 4 tests.
[       OK ] DeConvergenceMPI.Rastrigin (5 ms)
[----------] 4 tests from DeConvergenceMPI (28 ms total)

[----------] Global test environment tear-down
[==========] 4 tests from 1 test suite ran. (28 ms total)
[  PASSED  ] 4 tests.

```

