#ifndef GENETIC_TSP_THREADS_H
#define GENETIC_TSP_THREADS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <mutex>
#include <fstream>
#include <sstream>
#include <thread>
#include <functional>
#include "utimer.h"
#include "utils.h"

// Class that implements a map parallelization pattern
class ParallelMap {
private:
    int num_workers;

public:
    ParallelMap(int nw) : num_workers(nw) {}

    // Map over a vector of data
    template <typename T, typename Function, typename... Args>
    void execute(Function&& function, std::vector<T>& data, Args&&... args) {
        int chunk_size = data.size() / num_workers;

        std::vector<std::thread> threads;

        for (int i = 0; i < num_workers; i++) {
            // Calculate the chunks
            int start = i * chunk_size;
            int end;
            if (i == num_workers - 1) {
                end = data.size();
            } else {
                end = start + chunk_size;
            }
            threads.emplace_back([&function, &data, start, end, &args...]() {
                // Iterate over the chunk in each thread
                for (int j = start; j < end; j++) {
                    function(data[j], std::forward<Args>(args)...);
                }
            });
        }

        for (auto& thread : threads) {
            thread.join();
        }
    }

    // Map over a vector of input data and save the results in a vector of output data
    template <typename T, typename U, typename Function, typename... Args>
    void execute_save(Function&& function, std::vector<T>& input, std::vector<U>& output, Args&&... args) {
        int chunk_size = input.size() / num_workers;

        std::vector<std::thread> threads;

        for (int i = 0; i < num_workers; i++) {
            // Calculate the chunks
            int start = i * chunk_size;
            int end;
            if (i == num_workers - 1) {
                end = input.size();
            } else {
                end = start + chunk_size;
            }
            threads.emplace_back([&function, &input, &output, start, end, &args...]() {
                // Iterate over the chunk in each thread
                for (int j = start; j < end; j++) {
                    output[j] = function(input[j], std::forward<Args>(args)...);
                }
            });
        }

        for (auto& thread : threads) {
            thread.join();
        }
    }
};

// Function to execute the genetic algorithm using multiple threads
void experiment_threads(const int population_size, const int num_iterations, const float mutation_rate,
                        const float elitism_rate, const std::string citiesPth, const int num_workers);

#endif