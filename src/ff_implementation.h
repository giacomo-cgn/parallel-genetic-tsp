#ifndef GENETIC_TSP_FF_H
#define GENETIC_TSP_FF_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <mutex>
#include <fstream>
#include <sstream>
#include "utimer.h"
#include "utils.h"
#include "ff/parallel_for.hpp"

void experiment_ff(const int population_size, const int numIterations, const float mutation_rate,
                            const float elitism_rate, const std::string citiesPth, int numWorkers, bool recordInternalTimes,
                            bool printIterations);

#endif