#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <mutex>
#include <fstream>
#include <sstream>
#include "utimer.h"
#include "utils.h"
#include "ff_implementation.h"
#include <ff/parallel_for.hpp>



void experiment_ff(const int population_size, const int numIterations, const float mutation_rate,
                            const float elitism_rate, const std::string citiesPth, const int numWorkers, bool recordInternalTimes,
                            bool printIterations) {

    // Important variables
    std::vector<City> cities;
    std::vector<Chromosome> oldPopulation;
    std::vector<Chromosome> nextPopulation;
    std::vector<std::vector<float>> adjacencyMatrix;


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

    long totalTime;
    long distanceTime;
    long initializationTimeRandom;
    long initializationTimeEmpty;
    long evolutionTime = 0;
    // Internal times in evolutionTime
    long crossoverTime = 0;
    long mutationTime = 0;
    long fitnessTime = 0;

    // Initialize parallel_for
    ff::ParallelFor pf(numWorkers);


    std::vector<float> bestFitnesses(numIterations);

    {
        utimer timer(&totalTime);
    
    
        // Initialize adjacency matrix
        adjacencyMatrix.resize(cities.size());
        
        {
            utimer timer(&distanceTime);

            // Calculate the distance between each pair of cities and store it in an adjacency matrix  
            pf.parallel_for(0, cities.size(), 1, 0, [&](const long i) {   
                adjacencyMatrix[i] = generateDistanceRow(cities[i], cities);
            });
        }

        // Initialize random seed
        std::srand(std::time(nullptr));

        // Generate initial old population (random)
        oldPopulation.resize(population_size);
        
        {
            utimer timer(&initializationTimeRandom);

            ff::ParallelFor pf(numWorkers);
            pf.parallel_for(0, population_size, 1, 0, [&](const long i) {
                generateRandomChromosome(oldPopulation[i], cities, adjacencyMatrix);
            });
        }

        // Generate initial next population (empty)
        nextPopulation.resize(population_size);
        
        {
            utimer timer(&initializationTimeEmpty);

            pf.parallel_for(0, population_size, 1, 0, [&](const long i) {
                generateEmptyChromosome(nextPopulation[i], cities);
            });        
        }    

        // Start evolution iterations     
        for (int i = 0; i < numIterations; ++i) {
            int numBestParents = population_size * elitism_rate;
            // Sort the population in descending order based on fitness
            std::sort(oldPopulation.begin(), oldPopulation.end(), [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness;
            });

            if (printIterations) {
                std::cout << "Best fitness at iteration " << i-1 << ": " << oldPopulation[0].fitness << std::endl;
            }
            // Save best fitness
            bestFitnesses[i] = oldPopulation[0].fitness;

            long t;
            {
                utimer timer(&t);

                // Iterate over each chromosome in the next population and generate a child
                pf.parallel_for(0, population_size, 1, 0, [&](const long j) {
                    generateChild(nextPopulation[j], oldPopulation, numBestParents, mutation_rate, cities,
                                adjacencyMatrix, crossoverTime, mutationTime, fitnessTime, recordInternalTimes);
                });
            }
            evolutionTime += t;

            // Invert the populations for next iteration
            std::swap(oldPopulation, nextPopulation);
        }

    }

    std::cout << "######## FASTFLOW TIMES " << numWorkers <<" ########" << std::endl;

    
    std::cout << "DISTANCE MATRIX TIME: " << distanceTime << std::endl;
    std::cout << "POPULATION INITIALIZATION TIME RANDOM: " << initializationTimeRandom << std::endl;
    std::cout << "POPULATION INITIALIZATION TIME EMPTY: " << initializationTimeEmpty << std::endl;
    std::cout << "EVOLUTION TIME: " << evolutionTime << std::endl;
    std::cout << "TOTAL TIME: " << totalTime << std::endl;
    std::cout << "crossover time: " << crossoverTime << std::endl;
    std::cout << "mutation time: " << mutationTime << std::endl;
    std::cout << "fitness time: " << fitnessTime << std::endl;
    

    // Find final best fitness
    float maxFitness = oldPopulation[0].fitness;
    for (int i = 0; i < population_size; ++i) {
        if (oldPopulation[i].fitness > maxFitness) {
            maxFitness = oldPopulation[i].fitness;
        }
    }

    // Save the best fitness and times in a .csv in ../results
    std::ofstream resultsFile;
    resultsFile.open("../results/ff_results.csv", std::ios_base::app);
    // If not present, write header with column names
    if (resultsFile.tellp() == 0) {
        resultsFile << "numWorkers,maxFitness,totalTime,distanceTime,initializationTimeRandom,initializationTimeEmpty,evolutionTime,crossoverTime,mutationTime,fitnessTime" << std::endl;
    }
    resultsFile << numWorkers << "," << maxFitness << "," << totalTime << "," << distanceTime << "," << initializationTimeRandom << "," << initializationTimeEmpty << "," << evolutionTime << "," << crossoverTime << "," << mutationTime << "," << fitnessTime << std::endl;

    // Save the best fitnesses in a .csv in ../results. Each row has all the best fitnesses for a single run
    std::ofstream bestFitnessesFile;
    bestFitnessesFile.open("../results/best_fitnesses.csv", std::ios_base::app);
    bestFitnessesFile << "fastflow" << numWorkers << ",";
    for (int i = 0; i < numIterations; ++i) {
        bestFitnessesFile << bestFitnesses[i] << ",";
    }
    bestFitnessesFile << std::endl;
}
