#ifndef GENETIC_TSP_UTILS_H
#define GENETIC_TSP_UTILS_H

#include <vector>
#include <string>

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

    Chromosome();
};

// Function to calculate the distance between two cities
float calculateDistance(const City& city1, const City& city2);

// Function to generate one row of the adjacency matrix of distances between cities
std::vector<float> generateDistanceRow(const City& city, const std::vector<City>& cities);

// Function to get distance from the adjacency matrix
float getDistance(int cityIndex1, int cityIndex2, const std::vector<std::vector<float>>& adjacencyMatrix);

// Function to calculate the fitness of a chromosome
float calculateFitness(const Chromosome& chromosome, const std::vector<City>& cities, const std::vector<std::vector<float>>& adjacencyMatrix);

// Function to generate a random chromosome
void generateRandomChromosome(Chromosome& chromosome, const std::vector<City>& cities, const std::vector<std::vector<float>>& adjacencyMatrix);

// Generate empty chromosome
void generateEmptyChromosome(Chromosome& chromosome, const std::vector<City>& cities);

// Function to perform crossover between two parent chromosomes. It uses "Ordered Crossover (OX)"
void crossover(Chromosome& child, const Chromosome& parent1, const Chromosome& parent2, const std::vector<City>& cities);

// Mutate function: up to MUTATION_RATE*NUM_CITIES cities are swapped
void mutate(Chromosome& chromosome, const float mutation_rate, const std::vector<City>& cities);

// Function to generate a child chromosome from two parent chromosomes
void generateChild(Chromosome& child, const std::vector<Chromosome> oldPopulation, const int numBestParents,
                    const float mutation_rate, const std::vector<City>& cities, const std::vector<std::vector<float>>& adjacencyMatrix,
                    long* crossoverTime=nullptr, long* mutationTime=nullptr, long* fitnessTime=nullptr);

void printPopulation(const std::vector<Chromosome>& oldPopulation, std::string title, const std::vector<City>& cities);

#endif