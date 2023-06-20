#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <mutex>
#include <fstream>
#include <sstream>
#include "utimer.h"


// Constants
const int POPULATION_SIZE = 500;
const int NUM_ITERATIONS = 10;
const float MUTATION_RATE = 0.1;
const float ELITISM_RATE = 0.1;

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

// Global variables
std::vector<City> cities;
std::vector<Chromosome> oldPopulation;
std::vector<Chromosome> nextPopulation;

long crossoverTime = 0;
long mutationTime = 0;
long fitnessTime = 0;

// Function to calculate the distance between two cities
float calculateDistance(const City& city1, const City& city2) {
    float dx = city1.x - city2.x;
    float dy = city1.y - city2.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Function to generate the adjacency matrix of distances between cities
std::vector<std::vector<float>> generateAdjacencyMatrix(const std::vector<City>& cities) {
    std::vector<std::vector<float>> adjacencyMatrix(cities.size(), std::vector<float>(cities.size(), 0.0));
    for (int i = 0; i < cities.size(); ++i) {
        for (int j = 0; j < cities.size(); ++j) {
            adjacencyMatrix[i][j] = calculateDistance(cities[i], cities[j]);
        }
    }
    return adjacencyMatrix;
}

// Function to get distance from the adjacency matrix
float getDistance(const std::vector<std::vector<float>>& adjacencyMatrix, int cityIndex1, int cityIndex2) {
    return adjacencyMatrix[cityIndex1][cityIndex2];
}

// Function to calculate the fitness of a chromosome
float calculateFitness(const Chromosome& chromosome, const std::vector<std::vector<float>>& adjacencyMatrix) {
    float totalDistance = 0.0;
    for (int i = 0; i < cities.size() - 1; ++i) {
        int cityIndex1 = chromosome.path[i];
        int cityIndex2 = chromosome.path[i + 1];
        totalDistance += getDistance(adjacencyMatrix, cityIndex1, cityIndex2);
    }
    // Add the distance from the last city back to the first city
    totalDistance += getDistance(adjacencyMatrix, chromosome.path[cities.size() - 1], chromosome.path[0]);
    return 1.0 / totalDistance;
}

// Function to generate a random chromosome
Chromosome generateRandomChromosome(const std::vector<std::vector<float>>& adjacencyMatrix) {
    Chromosome chromosome;
    chromosome.path.resize(cities.size());
    for (int i = 0; i < cities.size(); ++i) {
        chromosome.path[i] = i;
    }
    std::random_shuffle(chromosome.path.begin() + 1, chromosome.path.end());
    chromosome.fitness = calculateFitness(chromosome, adjacencyMatrix);
    return chromosome;
}

// Function to perform crossover between two parent chromosomes. It uses "Ordered Crossover (OX)"
void crossover(Chromosome& child, const Chromosome& parent1, const Chromosome& parent2) {
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
void mutate(Chromosome& chromosome) {
   int numMutations = std::rand() % (int)(MUTATION_RATE * cities.size());
    for (int i = 0; i < numMutations; ++i) {
         int index1 = std::rand() % cities.size();
         int index2 = std::rand() % cities.size();
         std::swap(chromosome.path[index1], chromosome.path[index2]);
    }
}


// Function to generate a child chromosome from two parent chromosomes
void generateChild(Chromosome& child, const std::vector<Chromosome> oldPopulation, const int numBestParents, const std::vector<std::vector<float>>& adjacencyMatrix) {
    // Choose 2 random parents from oldPopulation up to numBestParents
    int parentIndex1 = std::rand() % numBestParents;
    int parentIndex2 = std::rand() % numBestParents;
    const Chromosome& parent1 = oldPopulation[parentIndex1];
    const Chromosome& parent2 = oldPopulation[parentIndex2];

    // Perform crossover and mutation to create a new child
    long t;
    {
        utimer timer(&t);
        crossover(child, parent1, parent2);

    }
    crossoverTime += t;
    {
        utimer timer(&t);
        mutate(child);
    }
    mutationTime += t;
    {
        utimer timer(&t);
        child.fitness = calculateFitness(child, adjacencyMatrix);
    }
    fitnessTime += t;
}
    
void printPopulation(const std::vector<Chromosome>& oldPopulation, std::string title) {
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

int main(int argc, char** argv) {

    // START INITIALIZATION

    // std::string citiesPth = "data/fi10639.tsp";
    // if (argc == 2 && std::string(argv[1]) == "--test") {
    //     citiesPth = "data/wi29.tsp";
    // }

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

    // Calculate the distance between each pair of cities and store it in an adjacency matrix
    std::vector<std::vector<float>> adjacencyMatrix = generateAdjacencyMatrix(cities);

    // Initialize random seed
    std::srand(std::time(nullptr));

    // Generate initial old population (random)
    long generationTime;
    {
        utimer timer(&generationTime);
        oldPopulation.resize(POPULATION_SIZE);
        for (int i = 0; i < POPULATION_SIZE; ++i) {
            oldPopulation[i] = generateRandomChromosome(adjacencyMatrix);
        }
        // Generate initial next population (empty)
        nextPopulation.resize(POPULATION_SIZE);
        for (int i = 0; i < POPULATION_SIZE; ++i) {
            Chromosome chromosome;
            chromosome.path.resize(cities.size());
            nextPopulation[i] = chromosome;
        }
    }    

    // Start evolution iterations
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
            generateChild(nextPopulation[j], oldPopulation, numBestParents, adjacencyMatrix);
        }

        // Invert the populations for next iteration
        std::swap(oldPopulation, nextPopulation);
    }
  
    std::cout << "GENERATION POPULATION TIME: " << generationTime << std::endl;
    std::cout << "CROSSOVER TIME: " << crossoverTime << std::endl;
    std::cout << "MUTATION TIME: " << mutationTime << std::endl;
    std::cout << "FITNESS TIME: " << fitnessTime << std::endl;

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
