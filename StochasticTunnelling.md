#  Stochastic Tunnelling Approach
Stochastic Tunnelling (ST) is an optimization algorithm designed for determining the global minima of complex and rugged energy landscapes. It is a generic physically motivated generalization of
simulated annealing. This approach circumvents the freezing problem which arises when the energy difference between “adjacent” local minima on the energy surface is much smaller than the energy of 
intervening transition states separating them.The physical idea behind the stochastic tunneling method is to allow the particle to “tunnel” forbidden regions, once it has been determined
that they are irrelevant for the low-energy properties of the problem. This can be accomplished by applying the transformation:

$$
f_{\text{STUN}}(x) = 1 - \exp[-\gamma(f(x) - f_0)],
$$

where $f_0$ is the lowest minimum encountered thus far. The effective potential preserves the locations of all minima, but maps the entire energy space from $f_0$ to the maximum
of the potential onto the interval [0, 1].




## Multithreading version (OMP)
### Main
> ⚙️ **Executable Overview** – This section explains how to run the main OpenMP Stochastic Tunnelling implementation.


#### Execution
An example invocation of the `main` executable for the OpenMP Stochastic Tunnelling implementation (all the detailed flag for main execution are in the main readme page):

```bash
./build/main -a stochastic_omp -d 2 -n 10 -i 1000 -f sphere
```

#### Sample output
```console

Minimum found:
  f(3.185893e-04, -3.450245e-03) = 1.529150e-07
  Total execution time: 0.024336 seconds
```

### Test
> ✅ **Test Suite Summary** – Automated tests  to validate correctness and convergence.


A comprehensive test suite is provided to verify the correctness and convergence properties of the Stochastic Tunnelling implementation. All tests are written in C++17 using the GoogleTest framework and exercise the algorithm on four classic benchmark functions.

#### Function tested 
| Test Name                         | Objective Function | Convergence Criterion                |
| --------------------------------- | ------------------ | ------------------------------------ |
| `TunnellingConvergence.Sphere`            | Sphere             | ‖f(x) – 0‖ ≤ 1 × 10⁻3                |
| `TunnellingConvergence.Rosenbrock`        | Rosenbrock         | ‖f(x) – 0‖  ≤ 1 × 10⁻2 |
| `TunnellingConvergence.Rastrigin`         | Rastrigin          | ‖f(x) – 0‖ ≤ 1 × 10⁻1               |
| `TunnellingConvergence.EuclideanDistance`         | Euclidean Distance         | ‖f(x) – 0‖ ≤ 1 × 10⁻2               |


#### Test Setup

- **Framework:** GoogleTest  
- **Dimensions:** 2  
- **Population size:** 100 candidates  
- **Max iterations:** 1 000  
- **Random seed:** 42  
- **Search bounds:** \[-10, 10\]  
- **Threads:** 5
- **Starting value of sigma (sigma_max)** = 1.0
- **Final value of sigma (sigma_min)** = 5.e-5
- **Gamma** = 0.0001
- **Beta_adjust_factor** = 0.9
- **Number of step used to compute avg of movements (tunnelling)** = 10
- **Threshold (beta_tresholding)** = 0.2
- **Frequency of best-position exchange among particles (time_step_updating)** = 100
  

Each test invokes:
```cpp
	const std::pair<std::vector<double>, double> result = 
		algorithm::run_multi_stochastic_tunnelling(
    dimensions,
    max_iterations,
    seed,
    lower_bound,
    upper_bound,
    sigma_max,
    sigma_min,
    s,
    gamma,
    beta_adjust_factor,
    false /*verbose*/,
    beta,
    tunnelling,
    beta_tresholding,
    num_positions,
    time_step_updating,
    num_threads);
```
#### Test result 
All four tests passed:

```console
[==========] Running 4 tests from 1 test suite.
[----------] Global test environment set-up.
[----------] 4 tests from TunnellingConvergence
[ RUN      ] TunnellingConvergence.Sphere
[       OK ] TunnellingConvergence.Sphere (53 ms)
[ RUN      ] TunnellingConvergence.Rosenbrock
[       OK ] TunnellingConvergence.Rosenbrock (24 ms)
[ RUN      ] TunnellingConvergence.Rastrigin
[       OK ] TunnellingConvergence.Rastrigin (28 ms)
[ RUN      ] TunnellingConvergence.EuclideanDistance
[       OK ] TunnellingConvergence.EuclideanDistance (42 ms)
[----------] 4 tests from TunnellingConvergence (149 ms total)

[----------] Global test environment tear-down
[==========] 4 tests from 1 test suite ran. (149 ms total)
[  PASSED  ] 4 tests.
```


#### Summary of test duration
| Test                   | Duration |
| ---------------------- | -------: |
| Sphere                 |    53 ms |
| Rosenbrock             |    24 ms |
| Rastrigin              |    28 ms |
| Euclidean Distance     |    42 ms
| **Total elapsed time** |   280 ms |



## Multiprocessing version (MPI)

> **Note:**
> Follow the instructions in the main readme page to compile properly to use MPI.
> This MPI-based implementation is provided for educational purposes. In practice, MPI shines on distributed systems or clusters, whereas on a single multi-core machine it can incur additional communication overhead. As a result, we do not include a full benchmark suite here (and the Google Benchmark integration currently causes conflicts). What remains is:
> 1. an MPI-enabled `main` that runs Differential Evolution across multiple processes, and  
> 2. the same convergence tests you’ve already seen (using GoogleTest).

### Test
You can run the tests like so:
```bash
 ./build/test_tunnelling_convergence_mpi
```

All the test are passed even in the multiprocessing version of the algorithm:
```console
[==========] Running 3 tests from 1 test suite.
[----------] Global test environment set-up.
[----------] 3 tests from GeneticConvergenceMPI
[ RUN      ] GeneticConvergenceMPI.Sphere

Minimum found:
  f(3.735023e-04, -2.134572e-03) = 4.695902e-06
  Total execution time: 0.040320 seconds

[       OK ] GeneticConvergenceMPI.Sphere (41 ms)
[ RUN      ] GeneticConvergenceMPI.Rosenbrock

Minimum found:
  f(9.923102e-01, 9.833503e-01) = 7.680052e-05
  Total execution time: 0.049206 seconds

[       OK ] GeneticConvergenceMPI.Rosenbrock (49 ms)
[ RUN      ] GeneticConvergenceMPI.Rastrigin

Minimum found:
  f(2.126418e-02, 7.994981e-03) = 1.022519e-01
  Total execution time: 0.087717 seconds

[       OK ] GeneticConvergenceMPI.Rastrigin (88 ms)
[----------] 3 tests from GeneticConvergenceMPI (180 ms total)

[----------] Global test environment tear-down
[==========] 3 tests from 1 test suite ran. (180 ms total)
[  PASSED  ] 3 tests.

```

