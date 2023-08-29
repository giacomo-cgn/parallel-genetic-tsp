#include <iostream>
#include "sequential.h"
#include "threads.h"
#include "ff_implementation.h"


int main(int, char**){

    // Constants
    const int POPULATION_SIZE = 1000;
    const int NUM_ITERATIONS = 50;
    const float MUTATION_RATE = 0.02;
    const float ELITISM_RATE = 0.1;
    const std::string CITIES_PTH = "../data/zi929.tsp";

    const std::vector<int> NUM_WORKER_LIST = {1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 
                                              24, 26, 28, 30, 32};

    std::cout << "######### STARTING EXPERIMENTS ##########\n";

    std::cout << "######### EXPERIMENT SEQUENTIAL ##########\n";
    experiment_sequential(POPULATION_SIZE, NUM_ITERATIONS, MUTATION_RATE, ELITISM_RATE, CITIES_PTH, false, false);

    for (int numWorkers : NUM_WORKER_LIST){
        std::cout << "######### EXPERIMENT THREADS:" << numWorkers << " ##########";
        experiment_threads(POPULATION_SIZE, NUM_ITERATIONS, MUTATION_RATE, ELITISM_RATE, CITIES_PTH, numWorkers, false, false);

        std::cout << "######### EXPERIMENT FASTFLOW:" << numWorkers << " ##########";
        experiment_ff(POPULATION_SIZE, NUM_ITERATIONS, MUTATION_RATE, ELITISM_RATE, CITIES_PTH, numWorkers, false, false);
    }

    return 0;
}
