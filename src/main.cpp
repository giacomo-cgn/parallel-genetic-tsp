#include <iostream>
#include "sequential.h"
#include "threads.h"
#include "ff_implementation.h"


int main(int, char**){

    // Constants
    const int POPULATION_SIZE = 1000;
    const int NUM_ITERATIONS = 100;
    const float MUTATION_RATE = 0.02;
    const float ELITISM_RATE = 0.1;
    const std::string CITIES_PTH = "../data/zi929.tsp";

    const std::vector<int> NUM_WORKER_LIST = {1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 
                                              24, 26, 28, 30, 32};

    std::cout << "######### STARTING EXPERIMENTS ##########\n";

    auto start1 = std::chrono::high_resolution_clock::now();
    std::cout << "######### EXPERIMENT SEQUENTIAL ##########\n";
    experiment_sequential(POPULATION_SIZE, NUM_ITERATIONS, MUTATION_RATE, ELITISM_RATE, CITIES_PTH, true, false);
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff1 = end1-start1;
    std::cout << "Total time seq: " << diff1.count() << " s\n";

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "######### EXPERIMENT SEQUENTIAL ##########\n";
    experiment_sequential(POPULATION_SIZE, NUM_ITERATIONS, MUTATION_RATE, ELITISM_RATE, CITIES_PTH, false, false);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end-start;
    std::cout << "Total time seq: " << diff.count() << " s\n";

    for (int numWorkers : NUM_WORKER_LIST){
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "######### EXPERIMENT THREADS:" << numWorkers << " ##########\n";
        experiment_threads(POPULATION_SIZE, NUM_ITERATIONS, MUTATION_RATE, ELITISM_RATE, CITIES_PTH, numWorkers, false, false);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end-start;
        std::cout << "Total time threads: " << diff.count() << " s\n";

        auto start_ff = std::chrono::high_resolution_clock::now();
        std::cout << "######### EXPERIMENT FASTFLOW:" << numWorkers << " ##########\n";
        experiment_ff(POPULATION_SIZE, NUM_ITERATIONS, MUTATION_RATE, ELITISM_RATE, CITIES_PTH, numWorkers, false, false);
        auto end_ff = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff_ff = end_ff-start_ff;
        std::cout << "Total time ff: " << diff_ff.count() << " s\n";
    }

    return 0;
}
