#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include "utimer.h"


// City struct
struct City {
    int id;
    float x;
    float y;
};

// Chromosome struct
struct Chromosome {
    std::vector<int> path;
    float fitness;

    Chromosome() : fitness(0.0) {}
};

// Function to calculate the distance between two cities
float calculateDistance(const City& city1, const City& city2) {
    float dx = city1.x - city2.x;
    float dy = city1.y - city2.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Function to generate one row of the adjacency matrix of distances between cities
std::vector<float> generateDistanceRow(const City& city, const std::vector<City>& cities) {
    std::vector<float> row(cities.size(), 0.0);
    for (int i = 0; i < cities.size(); ++i) {
        row[i] = calculateDistance(city, cities[i]);
    }
    return row;
}

// Function to get distance from the adjacency matrix
float getDistance(int cityIndex1, int cityIndex2, const std::vector<std::vector<float>>& adjacencyMatrix) {
    return adjacencyMatrix[cityIndex1][cityIndex2];
}

// Function to calculate the fitness of a chromosome
float calculateFitness(const Chromosome& chromosome, const std::vector<City>& cities, const std::vector<std::vector<float>>& adjacencyMatrix) {
    float totalDistance = 0.0;
    for (int i = 0; i < cities.size() - 1; ++i) {
        int cityIndex1 = chromosome.path[i];
        int cityIndex2 = chromosome.path[i + 1];
        totalDistance += getDistance(cityIndex1, cityIndex2, adjacencyMatrix);
    }
    // Add the distance from the last city back to the first city
    totalDistance += getDistance(chromosome.path[cities.size() - 1], chromosome.path[0], adjacencyMatrix);
    return 1.0 / totalDistance;
}

// Function to generate a random chromosome
void generateRandomChromosome(Chromosome& chromosome, const std::vector<City>& cities, const std::vector<std::vector<float>>& adjacencyMatrix) {
    chromosome.path.resize(cities.size());
    for (int i = 0; i < cities.size(); ++i) {
        chromosome.path[i] = i;
    }
    std::random_shuffle(chromosome.path.begin(), chromosome.path.end());
    chromosome.fitness = calculateFitness(chromosome, cities, adjacencyMatrix);
}

// Generate empty chromosome
void generateEmptyChromosome(Chromosome& chromosome, const std::vector<City>& cities) {
    chromosome.path.resize(cities.size());
    chromosome.fitness = 0.0;
}

// Function to perform crossover between two parent chromosomes. It uses "Ordered Crossover (OX)"
void crossover(Chromosome& child, const Chromosome& parent1, const Chromosome& parent2, const std::vector<City>& cities) {
    std::vector<bool> visited(cities.size(), false);

    // Select a random subset of the parent 1's path
    int randA = std::rand() % cities.size();
    int randB = std::rand() % cities.size();
    int startPos = std::min(randA, randB);
    int endPos = std::max(randA, randB);

    // Copy the subset of parent 1's path to the child's path
    for (int i = startPos; i <= endPos; ++i) {
        child.path[i] = parent1.path[i];
        visited[child.path[i]] = true;
    }

    // Fill the rest of the child's path with the remaining cities from parent 2, in order left to right
    int parent2Pos = 0;
    for (int i = 0; i < cities.size(); ++i) {
        if (parent2Pos == startPos)
            parent2Pos = endPos + 1;

        if (!visited[parent2.path[i]]) {
            child.path[parent2Pos] = parent2.path[i];
            visited[parent2.path[i]] = true;
            ++parent2Pos;
        }
    }
}

// Mutate function: up to MUTATION_RATE*NUM_CITIES cities are swapped
void mutate(Chromosome& chromosome, const float mutation_rate, const std::vector<City>& cities) {
   int numMutations = std::rand() % (int)(mutation_rate * cities.size());
    for (int i = 0; i < numMutations; ++i) {
         int index1 = std::rand() % cities.size();
         int index2 = std::rand() % cities.size();
         std::swap(chromosome.path[index1], chromosome.path[index2]);
    }
}


// Function to generate a child chromosome from two parent chromosomes
void generateChild(Chromosome& child, const std::vector<Chromosome> oldPopulation, const int numBestParents,
 const float mutation_rate, const std::vector<City>& cities, const std::vector<std::vector<float>>& adjacencyMatrix,
 long* crossoverTime=nullptr, long* mutationTime=nullptr, long* fitnessTime=nullptr) {
    // Choose 2 random parents from oldPopulation up to numBestParents
    int parentIndex1 = std::rand() % numBestParents;
    int parentIndex2 = std::rand() % numBestParents;
    const Chromosome& parent1 = oldPopulation[parentIndex1];
    const Chromosome& parent2 = oldPopulation[parentIndex2];


    if (crossoverTime != nullptr && mutationTime != nullptr && fitnessTime != nullptr) {
        // Perform crossover and mutation to create a new child recording times
        long t;
        {
            utimer timer(&t);
            crossover(child, parent1, parent2, cities);

        }
        crossoverTime += t;
        {
            utimer timer(&t);
            mutate(child, mutation_rate, cities);
        }
        mutationTime += t;
        {
            utimer timer(&t);
            child.fitness = calculateFitness(child, cities, adjacencyMatrix);
        }
        fitnessTime += t;
    }
    else {
        // Perform crossover and mutation to create a new child without recording times
        crossover(child, parent1, parent2, cities);
        mutate(child, mutation_rate, cities);
        child.fitness = calculateFitness(child, cities, adjacencyMatrix);
    }
}
    
void printPopulation(const std::vector<Chromosome>& oldPopulation, std::string title, const std::vector<City>& cities) {
    std::cout << title << std::endl;

    // Print all the paths next to their fitness
    for (const auto& chromosome : oldPopulation) {
        std::cout << "Path: ";
        for (int i = 0; i < cities.size(); ++i) {
            std::cout << chromosome.path[i] << " ";
        }
        std::cout << "Fitness: " << chromosome.fitness << std::endl;
    }
}