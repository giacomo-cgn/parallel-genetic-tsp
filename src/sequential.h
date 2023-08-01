#ifndef GENETIC_TSP_SEQUENTIAL_H
#define GENETIC_TSP_SEQUENTIAL_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <mutex>
#include <fstream>
#include <sstream>
#include "utimer.h"
#include "utils.h"
#include "sequential.h"

const int POPULATION_SIZE = 1000;
const int NUM_ITERATIONS = 50;
const float MUTATION_RATE = 0.02;
const float ELITISM_RATE = 0.1;

void experiment_sequential(const int population_size, const int num_iterations, const float mutation_rate,
                            const float elitism_rate, const std::string citiesPth);

#endif