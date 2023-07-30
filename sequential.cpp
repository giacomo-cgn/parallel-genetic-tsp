#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <mutex>
#include <fstream>
#include <sstream>
#include "utimer.h"
#include "utils.h"


// Constants
const int POPULATION_SIZE = 1000;
const int NUM_ITERATIONS = 50;
const float MUTATION_RATE = 0.02;
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
 
    // Initialize adjacency matrix
    long distanceTime;
    {
        utimer timer(&distanceTime);
        adjacencyMatrix.resize(cities.size());
        // Calculate the distance between each pair of cities and store it in an adjacency matrix    
        for (int i = 0; i < cities.size(); ++i) {     
            adjacencyMatrix[i] = generateDistanceRow(cities[i], cities);
        }
    }

    // Initialize random seed
    std::srand(std::time(nullptr));

    // Generate initial old population (random)
    long initializationTimeRandom;
    {
        utimer timer(&initializationTimeRandom);
        oldPopulation.resize(POPULATION_SIZE);
        for (int i = 0; i < POPULATION_SIZE; ++i) {
            generateRandomChromosome(oldPopulation[i], cities, adjacencyMatrix);
        }
    }

    // Generate initial next population (empty)
    long initializationTimeEmpty;
    {
        utimer timer(&initializationTimeEmpty);
        nextPopulation.resize(POPULATION_SIZE);
        for (int i = 0; i < POPULATION_SIZE; ++i) {
            generateEmptyChromosome(nextPopulation[i], cities);
        }
        
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
            for (int j = 0; j < POPULATION_SIZE; ++j) {
                generateChild(nextPopulation[j], oldPopulation, numBestParents, MUTATION_RATE, cities,
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

    // Print the best fitness and path at the end in the old population
    std::cout << "Best fitness at final iteration: " << oldPopulation[0].fitness << std::endl;
    std::cout << "Best path: ";
    for (int i = 0; i < cities.size(); ++i) {
        std::cout << oldPopulation[0].path[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
