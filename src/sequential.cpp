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


void experiment_sequential(const int population_size, const int num_iterations, const float mutation_rate,
                            const float elitism_rate, const std::string citiesPth) {

    // Important variables
    std::vector<City> cities;
    std::vector<Chromosome> oldPopulation;
    std::vector<Chromosome> nextPopulation;
    std::vector<std::vector<float>> adjacencyMatrix;

    long crossoverTime = 0;
    long mutationTime = 0;
    long fitnessTime = 0;

    // START INITIALIZATION

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
 
    adjacencyMatrix.resize(cities.size());
    // Initialize adjacency matrix
    long distanceTime;
    {
        utimer timer(&distanceTime);
        // Calculate the distance between each pair of cities and store it in an adjacency matrix    
        for (int i = 0; i < cities.size(); ++i) {     
            adjacencyMatrix[i] = generateDistanceRow(cities[i], cities);
        }
    }

    // Initialize random seed
    std::srand(std::time(nullptr));

    // Generate initial old population (random)
    oldPopulation.resize(population_size);
    long initializationTimeRandom;
    {
        utimer timer(&initializationTimeRandom);
        for (int i = 0; i < population_size; ++i) {
            generateRandomChromosome(oldPopulation[i], cities, adjacencyMatrix);
        }
    }

    // Generate initial next population (empty)
    nextPopulation.resize(population_size);
    long initializationTimeEmpty;
    {
        utimer timer(&initializationTimeEmpty);
        for (int i = 0; i < population_size; ++i) {
            generateEmptyChromosome(nextPopulation[i], cities);
        }
        
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
            for (int j = 0; j < population_size; ++j) {
                generateChild(nextPopulation[j], oldPopulation, numBestParents, mutation_rate, cities,
                            adjacencyMatrix, &crossoverTime, &mutationTime, &fitnessTime);
            }

            // Invert the populations for next iteration
            std::swap(oldPopulation, nextPopulation);
        }
    }

    std::cout << "######## SEQUENTIAL TIMES ########" << std::endl;

    
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

    // // Print the best fitness and path at the end in the old population
    // std::cout << "Best fitness at final iteration: " << oldPopulation[0].fitness << std::endl;
    // std::cout << "Best path: ";
    // for (int i = 0; i < cities.size(); ++i) {
    //     std::cout << oldPopulation[0].path[i] << " ";
    // }
    // std::cout << std::endl;

}
