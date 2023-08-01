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
#include "threads.h"


void experiment_threads(const int population_size, const int num_iterations, const float mutation_rate,
                            const float elitism_rate, const std::string citiesPth, const int num_workers) {

    // Important variables
    std::vector<City> cities;
    std::vector<Chromosome> oldPopulation;
    std::vector<Chromosome> nextPopulation;
    std::vector<std::vector<float>> adjacencyMatrix;

    long crossoverTime = 0;
    long mutationTime = 0;
    long fitnessTime = 0;

    // START INITIALIZATION

    // Initialize parallel map
    ParallelMap parMap(num_workers);

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
        oldPopulation.resize(population_size);
        parMap.execute(generateRandomChromosome, oldPopulation, cities, adjacencyMatrix);
        
    }

    // Generate initial next population (empty)
    long initializationTimeEmpty;
    {
        utimer timer(&initializationTimeEmpty);
        nextPopulation.resize(population_size);

        parMap.execute(generateEmptyChromosome, nextPopulation, cities);
    }    

    // Start evolution iterations
    long evolutionTime;
    {
        utimer timer(&evolutionTime);

        for (int i = 0; i < num_iterations; ++i) {
            int numBestParents = population_size * elitism_rate;
            // Sort the population in descending order based on fitness
            std::sort(oldPopulation.begin(), oldPopulation.end(), [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness;
            });

            // print the best fitness in the old population
            std::cout << "Best fitness at iteration " << i-1 << ": " << oldPopulation[0].fitness << std::endl;

            // Iterate over each chromosome in the next population and generate a child
            parMap.execute(generateChild, nextPopulation, oldPopulation, numBestParents, mutation_rate, cities,
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
    // std::cout << "Best fitness at final iteration: " << oldPopulation[0].fitness << std::endl;
    // std::cout << "Best path: ";
    // for (int i = 0; i < cities.size(); ++i) {
    //     std::cout << oldPopulation[0].path[i] << " ";
    // }
    // std::cout << std::endl;

}