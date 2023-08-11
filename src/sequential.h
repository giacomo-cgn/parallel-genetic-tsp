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

void experiment_sequential(const int population_size, const int num_iterations, const float mutation_rate,
                            const float elitism_rate, const std::string citiesPth, bool recordInternalTimes);

#endif