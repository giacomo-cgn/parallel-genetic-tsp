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
        ParallelMap(int nw)
            : num_workers(nw) {}

        // Map over a vector of data
        template <typename T, typename Function, typename... Args>
        void execute(Function&& function, std::vector<T>& data, Args... args) {
            int chunk_size = data.size() / num_workers;

            std::vector<std::thread> threads;

            for (int i = 0; i < num_workers; i++) {
                // Calculate the chunks
                int start = i * chunk_size;
                int end;
                if (i == num_workers - 1) {
                    end = data.size();
                }
                else {
                    end = start + chunk_size;
                }
                threads.emplace_back([&function, &data, start, end, args...]() {
                    // Iterate over the chunk in each thread
                    for (int j = start; j < end; j++) {
                        function(data[j], args...);
                    }
                });
            }

            for (auto& thread : threads) {
                thread.join();
            }
        }

        // Map over a vector of input data and and save the results in a vector of output data
        template <typename T, typename U, typename Function, typename... Args>
        void execute_save(Function&& function, std::vector<T>& input, std::vector<U>& output, Args... args) {
            int chunk_size = input.size() / num_workers;

            std::vector<std::thread> threads;

            for (int i = 0; i < num_workers; i++) {
                // Calculate the chunks
                int start = i * chunk_size;
                int end;
                if (i == num_workers - 1) {
                    end = input.size();
                }
                else {
                    end = start + chunk_size;
                }
                threads.emplace_back([&function, &input, &output, start, end, args...]() {
                    // Iterate over the chunk in each thread
                    for (int j = start; j < end; j++) {
                        output[j] = function(input[j], args...);
                    }
                });
            }

            for (auto& thread : threads) {
                thread.join();
            }
        }
};

// Constants
const int POPULATION_SIZE = 500;
const int NUM_ITERATIONS = 100;
const float MUTATION_RATE = 0.1;
const float ELITISM_RATE = 0.1;



int main(int argc, char** argv) {

    // Main variables
    std::vector<City> cities;
    std::vector<Chromosome> oldPopulation;
    std::vector<Chromosome> nextPopulation;
    std::vector<std::vector<float>> adjacencyMatrix;

    long crossoverTime = 0;
    long mutationTime = 0;
    long fitnessTime = 0;

    // START INITIALIZATION

    int num_workers = 8;
    // If -nw flag is passed, use the next argument as the number of workers
    if (argc > 1 && std::string(argv[1]) == "-nw") {
        num_workers = std::stoi(argv[2]);
    }
    // Initialize parallel map
    ParallelMap parMap(num_workers);

    // Path to cities data file
    std::string citiesPth = "data/zi929.tsp";

    // Read cities from file. Each row contains the x and y coordinates of a city separated by a space.
    std::ifstream file(citiesPth);
    std::string line;
    for (int i = 0; i < 7; ++i) {
        std::getline(file, line);
    }
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int id;
        float x, y;
        if (!(iss >> id >> x >> y)) {
            break;
        }
        City city;
        city.id = id - 1;
        city.x = x;
        city.y = y;
        cities.push_back(city);
    }

    long distanceTime;
    {
        utimer timer(&distanceTime);
        // Initialize adjacency matrix
        adjacencyMatrix.resize(cities.size());
        // Calculate the distance between each pair of cities and store it in an adjacency matrix    
        parMap.execute_save(generateDistanceRow, cities, adjacencyMatrix, cities);
    }


    // Initialize random seed
    std::srand(std::time(nullptr));

    // Generate initial old population (random)
    long initializationTimeRandom;
    {   
        utimer timer(&initializationTimeRandom);
        oldPopulation.resize(POPULATION_SIZE);
        parMap.execute(generateRandomChromosome, oldPopulation, cities, adjacencyMatrix);
        
    }

    // Generate initial next population (empty)
    long initializationTimeEmpty;
    {
        utimer timer(&initializationTimeEmpty);
        nextPopulation.resize(POPULATION_SIZE);

        parMap.execute(generateEmptyChromosome, nextPopulation, cities);
    }    

    // Start evolution iterations
    long evolutionTime;
    {
        utimer timer(&evolutionTime);

        for (int i = 0; i < NUM_ITERATIONS; ++i) {
            int numBestParents = POPULATION_SIZE * ELITISM_RATE;
            // Sort the population in descending order based on fitness
            std::sort(oldPopulation.begin(), oldPopulation.end(), [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness;
            });

            // print the best fitness in the old population
            std::cout << "Best fitness at iteration " << i-1 << ": " << oldPopulation[0].fitness << std::endl;

            // Iterate over each chromosome in the next population and generate a child
            parMap.execute(generateChild, nextPopulation, oldPopulation, numBestParents, MUTATION_RATE, cities,
                        adjacencyMatrix, &crossoverTime, &mutationTime, &fitnessTime);
    

            // Invert the populations for next iteration
            std::swap(oldPopulation, nextPopulation);
        }
    }

    std::cout << "######## THREAD TIMES ########" << std::endl;
  
    std::cout << "DISTANCE MATRIX TIME: " << distanceTime << std::endl;
    std::cout << "POPULATION INITIALIZATION TIME RANDOM: " << initializationTimeRandom << std::endl;
    std::cout << "POPULATION INITIALIZATION TIME EMPTY: " << initializationTimeEmpty << std::endl;
    std::cout << "EVOLUTION TIME: " << evolutionTime << std::endl;

    std::cout << "crossover time: " << crossoverTime << std::endl;
    std::cout << "mutation time: " << mutationTime << std::endl;
    std::cout << "fitness time: " << fitnessTime << std::endl;

    // Final sort of the population
    std::sort(oldPopulation.begin(), oldPopulation.end(), [](const Chromosome& a, const Chromosome& b) {
        return a.fitness > b.fitness;
    });

    // Print the best fitness and path at the end in the old population
    std::cout << "Best fitness at final iteration: " << oldPopulation[0].fitness << std::endl;
    std::cout << "Best path: ";
    for (int i = 0; i < cities.size(); ++i) {
        std::cout << oldPopulation[0].path[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
