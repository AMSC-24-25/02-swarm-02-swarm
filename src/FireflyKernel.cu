#ifdef ENABLE_CUDA


#include <curand_kernel.h>
#include <cmath>
#include "FireflyAlgorithm_Cuda.h"  // Declaration for host-callable wrapper

// Basic objective functions used for testing
// is kinda limit of CUDA (everytime you need a new function, you need to add here, for now)
__device__ double evalObjective(
    const double* x, int dimensions, int type
) {
    if (type == 0) { // Sphere function
        double sum = 0.0;
        for (int d = 0; d < dimensions; ++d)
            sum += x[d] * x[d];
        return sum;
    } else if (type == 1) { // Rastrigin function
        double A = 10.0, res = A * dimensions;
        for (int d = 0; d < dimensions; ++d)
            res += x[d]*x[d] - A * cos(2 * M_PI * x[d]);
        return res;
    } else if (type == 2) { // Rosenbrock function
        double sum = 0.0;
        for (int d = 0; d < dimensions - 1; ++d) {
            double t1 = (x[d + 1] - x[d] * x[d]);
            double t2 = (1 - x[d]);
            sum += 100 * t1 * t1 + t2 * t2;
        }
        return sum;
    }
    return 0.0; // fallback
}

// Kernel: evaluate objective function for each firefly
__global__ void evaluateBrightnessCUDA(
    double* positions, double* brightness,
    int numFireflies, int dimensions, int objectiveType
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numFireflies) return;
    brightness[i] = evalObjective(&positions[i * dimensions], dimensions, objectiveType);
}

// Kernel: update firefly positions using shared memory
__global__ void updateFirefliesCUDA_shared(
    double* positions, double* brightness,
    int numFireflies, int dimensions,
    double alpha, double beta, double gamma,
    unsigned int seed,
    double* new_positions
) {
    // Shared memory layout: [positions | brightness]
    extern __shared__ double shared[];
    double* shared_positions = shared;
    double* shared_brightness = (double*)&shared_positions[blockDim.x * dimensions];

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numFireflies) return;

    // Setup per-thread RNG
    curandState state;
    curand_init(seed, i, 0, &state);

    // Init: copy current position to output buffer
    for (int d = 0; d < dimensions; ++d)
        new_positions[i * dimensions + d] = positions[i * dimensions + d];

    // Iterate over all blocks (tiling)
    for (int block_j = 0; block_j * blockDim.x < numFireflies; ++block_j) {
        int j_global = block_j * blockDim.x + threadIdx.x;

        // Load batch into shared memory
        if (j_global < numFireflies) {
            for (int d = 0; d < dimensions; ++d)
                shared_positions[threadIdx.x * dimensions + d] = positions[j_global * dimensions + d];
            shared_brightness[threadIdx.x] = brightness[j_global];
        } else {
            // Fill with dummy data
            for (int d = 0; d < dimensions; ++d)
                shared_positions[threadIdx.x * dimensions + d] = 0.0;
            shared_brightness[threadIdx.x] = 1e20; // high dummy value
        }

        __syncthreads();

        // Update logic for firefly i
        for (int j_local = 0; j_local < blockDim.x; ++j_local) {
            int j = block_j * blockDim.x + j_local;
            if (j >= numFireflies) continue;

            if (shared_brightness[j_local] < brightness[i]) {
                // Compute Euclidean distance r
                double r = 0.0;
                for (int d = 0; d < dimensions; ++d) {
                    double diff = positions[i * dimensions + d] - shared_positions[j_local * dimensions + d];
                    r += diff * diff;
                }
                r = sqrt(r);
                double beta_eff = beta * exp(-gamma * r * r);

                // Move firefly i towards j with attractiveness and noise
                for (int d = 0; d < dimensions; ++d) {
                    double randNoise = 2.0 * curand_uniform_double(&state) - 1.0;
                    double xi = positions[i * dimensions + d];
                    double xj = shared_positions[j_local * dimensions + d];
                    new_positions[i * dimensions + d] += beta_eff * (xj - xi) + alpha * randNoise;
                }
            }
        }

        __syncthreads();
    }
}

// Host-callable wrapper to evaluate brightness
extern "C" void launchEvaluateBrightnessCUDA(
    double* d_positions, double* d_brightness,
    int numFireflies, int dimensions, int objectiveType,
    int threadsPerBlock
) {
    int blocksPerGrid = (numFireflies + threadsPerBlock - 1) / threadsPerBlock;
    evaluateBrightnessCUDA<<<blocksPerGrid, threadsPerBlock>>>(
        d_positions, d_brightness, numFireflies, dimensions, objectiveType
    );
    cudaDeviceSynchronize();
}

// Host-callable wrapper to update fireflies
extern "C" void launchUpdateFirefliesCUDA(
    double* d_positions, double* d_brightness,
    int numFireflies, int dimensions,
    double alpha, double beta, double gamma,
    unsigned int seed,
    int threadsPerBlock
) {
    int blocksPerGrid = (numFireflies + threadsPerBlock - 1) / threadsPerBlock;
    size_t sharedMemSize = threadsPerBlock * (dimensions * sizeof(double) + sizeof(double));

    // Allocate temporary buffer
    double* d_new_positions;
    cudaMalloc(&d_new_positions, numFireflies * dimensions * sizeof(double));

    updateFirefliesCUDA_shared<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(
        d_positions, d_brightness,
        numFireflies, dimensions,
        alpha, beta, gamma,
        seed,
        d_new_positions
    );

    // Copy back updated positions
    cudaMemcpy(d_positions, d_new_positions, numFireflies * dimensions * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaFree(d_new_positions);
    cudaDeviceSynchronize();
}

#endif // ENABLE_CUDA
