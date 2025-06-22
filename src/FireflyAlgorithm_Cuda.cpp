//#ifdef ENABLE_CUDA

#include "FireflyAlgorithm_Cuda.h"
#include <cuda_runtime.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include <iomanip>

#include <limits>
#include <omp.h>

// Constructor: calls base constructor (which already calls initializeFireflies)
FireflyAlgorithm_Cuda::FireflyAlgorithm_Cuda(int nFireflies, int dim, double a, double b, double g,double lower, double upper, size_t s)
    : FireflyAlgorithm(nFireflies,dim, a, b, g, lower, upper, s) {
    // nothing else needed here
}

// Kernel wrappers (defined in .cu file)
extern "C" void launchUpdateFirefliesCUDA(
    double* d_positions, double* d_brightness,
    int numFireflies, int dimensions,
    double alpha, double beta, double gamma,
    unsigned int seed,
    double lower_bound, double upper_bound,
    int threadsPerBlock
);

extern "C" void launchEvaluateBrightnessCUDA(
    double* d_positions, double* d_brightness,
    int numFireflies, int dimensions, int objectiveType,
    int threadsPerBlock
);

// Set objective function type (Sphere, Rastrigin, etc.)
void FireflyAlgorithm_Cuda::setObjectiveFunction(ObjectiveType type) {
    this->objectiveType = type;
}

void FireflyAlgorithm_Cuda::setObjectiveFunction(std::function<double(const std::vector<double>&)> func) {
    objectiveFunctionCPU = func;
    fitnessMode = FITNESS_CPU;
}


// Main optimization loop using CUDA
std::vector<double> FireflyAlgorithm_Cuda::optimize(int maxIterations) {
	auto start = std::chrono::high_resolution_clock::now();

    const int sizePos = numFireflies * dimensions;
    const size_t sizePosBytes = sizePos * sizeof(double);
    const size_t sizeBrightBytes = numFireflies * sizeof(double);

    std::vector<double> h_positions(sizePos);
    std::vector<double> h_brightness(numFireflies);

    // Flatten initial positions
    for (int i = 0; i < numFireflies; ++i) {
        const std::vector<double>& pos = fireflies[i].getPosition();
        for (int d = 0; d < dimensions; ++d)
            h_positions[i * dimensions + d] = pos[d];
    }

    // Allocate GPU memory
    double* d_positions;
    double* d_brightness;
    cudaMalloc(&d_positions, sizePosBytes);
    cudaMalloc(&d_brightness, sizeBrightBytes);
    cudaMemcpy(d_positions, h_positions.data(), sizePosBytes, cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;

    for (int iter = 0; iter < maxIterations; ++iter) {

        if (fitnessMode == FITNESS_CPU) {
            // 1. Copia le posizioni dalla GPU all’host
            cudaMemcpy(h_positions.data(), d_positions, sizePosBytes, cudaMemcpyDeviceToHost);

            // 2. Valuta la funzione obiettivo su CPU
#pragma omp parallel for
            for (int i = 0; i < numFireflies; ++i) {
                std::vector<double> fireflyPos(dimensions);
                for (int d = 0; d < dimensions; ++d)
                    fireflyPos[d] = h_positions[i * dimensions + d];
                h_brightness[i] = objectiveFunctionCPU(fireflyPos);
            }

            // 3. Copia i valori di brightness sulla GPU
            cudaMemcpy(d_brightness, h_brightness.data(), sizeBrightBytes, cudaMemcpyHostToDevice);

            // 4. Aggiorna le posizioni sulla GPU
            launchUpdateFirefliesCUDA(
                d_positions, d_brightness,
                numFireflies, dimensions,
                alpha, beta, gamma,
                static_cast<unsigned int>(time(NULL)) + iter,
                lower_bound, upper_bound,
                threadsPerBlock
            );

        	// Trova bestIndex
        	int bestIndex = 0;
        	for (int i = 1; i < numFireflies; ++i) {
        		if (h_brightness[i] < h_brightness[bestIndex])
        			bestIndex = i;
        	}

        	//  Stampa risultato dell’iterazione
        	std::cout << "Iteration n. " << iter + 1 << " / " << maxIterations << "\n";
        	std::cout << "  Current minimum:\n";
        	std::cout << "  f(";
        	for (int d = 0; d < dimensions; ++d) {
        		std::cout << std::scientific << std::setprecision(6)
						  << h_positions[bestIndex * dimensions + d];
        		if (d < dimensions - 1) std::cout << ", ";
        	}
        	std::cout << ") = " << std::scientific << std::setprecision(6)
					  << h_brightness[bestIndex] << std::endl;


        } else {
            // Caso GPU puro: come ora
            launchEvaluateBrightnessCUDA(
                d_positions, d_brightness,
                numFireflies, dimensions,
                objectiveType,
                threadsPerBlock
            );

            launchUpdateFirefliesCUDA(
                d_positions, d_brightness,
                numFireflies, dimensions,
                alpha, beta, gamma,
                static_cast<unsigned int>(time(NULL)) + iter,
                lower_bound, upper_bound,
                threadsPerBlock
            );
        }
    }

    // Resto invariato: copia finale e scelta del best
    cudaMemcpy(h_positions.data(), d_positions, sizePosBytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_brightness.data(), d_brightness, sizeBrightBytes, cudaMemcpyDeviceToHost);

    cudaFree(d_positions);
    cudaFree(d_brightness);

    for (int i = 0; i < numFireflies; ++i) {
        std::vector<double> pos(dimensions);
        for (int d = 0; d < dimensions; ++d)
            pos[d] = h_positions[i * dimensions + d];
        fireflies[i].setPosition(pos);
        fireflies[i].setBrightness(h_brightness[i]);
    }

    int bestIndex = 0;
    for (int i = 1; i < numFireflies; ++i) {
	    if (fireflies[i].getBrightness() < fireflies[bestIndex].getBrightness())
	    	bestIndex = i;
    }


	auto end = std::chrono::high_resolution_clock::now();
	double elapsed = std::chrono::duration<double>(end - start).count();
	std::cout << "Total execution time: " << std::fixed << std::setprecision(6)
			  << elapsed << " seconds" << std::endl;
	std::cout << std::endl;

    return fireflies[bestIndex].getPosition();
}

//#endif // ENABLE_CUDA