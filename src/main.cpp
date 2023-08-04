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
    const int NUM_WORKERS = 8;

    std::cout << "######### STARTING EXPERIMENTS ##########\n";

    std::cout << "######### EXPERIMENT SEQUENTIAL ##########\n";
    experiment_sequential(POPULATION_SIZE, NUM_ITERATIONS, MUTATION_RATE, ELITISM_RATE, CITIES_PTH);

    std::cout << "######### EXPERIMENT THREADS ##########\n";
    experiment_threads(POPULATION_SIZE, NUM_ITERATIONS, MUTATION_RATE, ELITISM_RATE, CITIES_PTH, NUM_WORKERS);

    std::cout << "######### EXPERIMENT FASTFLOW ##########\n";
    experiment_ff(POPULATION_SIZE, NUM_ITERATIONS, MUTATION_RATE, ELITISM_RATE, CITIES_PTH, NUM_WORKERS);

    return 0;
}
