#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <thread>
#include <mutex>

// Constants
const int NUM_CITIES = 10;
const int POPULATION_SIZE = 100;
const int NUM_ITERATIONS = 20;
const int MAX_COORDINATE = 1000;
const float MUTATION_RATE = 0.2;
const float ELITISM_RATE = 0.2;

// City struct
struct City {
    int id;
    int x;
    int y;
};

// Chromosome struct
struct Chromosome {
    std::vector<int> path;
    double fitness;

    Chromosome() : fitness(0.0) {}
};

// Global variables
std::vector<City> cities;
std::vector<Chromosome> oldPopulation;
std::vector<Chromosome> nextPopulation;    

// Function to calculate the distance between two cities
double calculateDistance(const City& city1, const City& city2) {
    double dx = city1.x - city2.x;
    double dy = city1.y - city2.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Function to calculate the fitness of a chromosome
double calculateFitness(const Chromosome& chromosome) {
    double totalDistance = 0.0;
    for (int i = 0; i < NUM_CITIES - 1; ++i) {
        int cityIndex1 = chromosome.path[i];
        int cityIndex2 = chromosome.path[i + 1];
        totalDistance += calculateDistance(cities[cityIndex1], cities[cityIndex2]);
    }
    // Add the distance from the last city back to the first city
    totalDistance += calculateDistance(cities[chromosome.path[NUM_CITIES - 1]], cities[chromosome.path[0]]);
    return 1.0 / totalDistance;
}

// Function to generate a random chromosome
Chromosome generateRandomChromosome() {
    Chromosome chromosome;
    chromosome.path.resize(NUM_CITIES);
    for (int i = 0; i < NUM_CITIES; ++i) {
        chromosome.path[i] = i;
    }
    std::random_shuffle(chromosome.path.begin() + 1, chromosome.path.end());
    chromosome.fitness = calculateFitness(chromosome);
    return chromosome;
}

// Function to perform crossover between two parent chromosomes
void crossover(Chromosome& child, const Chromosome& parent1, const Chromosome& parent2) {
    std::vector<bool> visited(NUM_CITIES, false);

    // Select a random subset of the parent 1's path
    int randA = std::rand() % NUM_CITIES;
    int randB = std::rand() % NUM_CITIES;
    int startPos = std::min(randA, randB);
    int endPos = std::max(randA, randB);

    // Copy the subset of parent 1's path to the child's path
    for (int i = startPos; i <= endPos; ++i) {
        child.path[i] = parent1.path[i];
        visited[child.path[i]] = true;
    }

    // Fill the rest of the child's path with the remaining cities from parent 2
    int parent2Pos = 0;
    for (int i = 0; i < NUM_CITIES; ++i) {
        if (parent2Pos == startPos)
            parent2Pos = endPos + 1;

        if (!visited[parent2.path[i]]) {
            child.path[parent2Pos] = parent2.path[i];
            visited[parent2.path[i]] = true;
            ++parent2Pos;
        }
    }
}


// Mutate function: up to MUTATION_RATE cities are swapped
void mutate(Chromosome& chromosome) {
   int numMutations = std::rand() % (int)(MUTATION_RATE * NUM_CITIES);
    for (int i = 0; i < numMutations; ++i) {
         int index1 = std::rand() % NUM_CITIES;
         int index2 = std::rand() % NUM_CITIES;
         std::swap(chromosome.path[index1], chromosome.path[index2]);
    }
}


// Function to generate a child chromosome from two parent chromosomes
void generateChild(Chromosome& child, const std::vector<Chromosome> oldPopulation, const int numBestParents) {
    // Choose 2 random parents from oldPopulation up to numBestParents
    int parentIndex1 = std::rand() % numBestParents;
    int parentIndex2 = std::rand() % numBestParents;
    const Chromosome& parent1 = oldPopulation[parentIndex1];
    const Chromosome& parent2 = oldPopulation[parentIndex2];

    // Perform crossover and mutation to create a new child
    crossover(child, parent1, parent2);
    mutate(child);
    child.fitness = calculateFitness(child);

}
    
void printPopulation(const std::vector<Chromosome>& oldPopulation, std::string title) {
    std::cout << title << std::endl;

    // Print all the paths next to their fitness
    for (const auto& chromosome : oldPopulation) {
        std::cout << "Path: ";
        for (int i = 0; i < NUM_CITIES; ++i) {
            std::cout << chromosome.path[i] << " ";
        }
        std::cout << "Fitness: " << chromosome.fitness << std::endl;
    }
}

int main() {
    // Initialize random seed
    std::srand(std::time(nullptr));

    // Generate random cities. x and y coordinates are ints between 0 and MAX_COORDINATE
    for (int i = 0; i < NUM_CITIES; ++i) {
        City city;
        city.id = i;
        city.x = std::rand() % MAX_COORDINATE;
        city.y = std::rand() % MAX_COORDINATE;
        cities.push_back(city);
    }

    // Generate initial old population
    oldPopulation.resize(POPULATION_SIZE);
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        oldPopulation[i] = generateRandomChromosome();
    }
    // Generate initial next population (empty)
    nextPopulation.resize(POPULATION_SIZE);
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        Chromosome chromosome;
        chromosome.path.resize(NUM_CITIES);
        nextPopulation[i] = chromosome;
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
            generateChild(nextPopulation[j], oldPopulation, numBestParents);
        }

        // Invert the populations for next iteration
        std::swap(oldPopulation, nextPopulation);
    }
  

    // Final sort of the population
    std::sort(oldPopulation.begin(), oldPopulation.end(), [](const Chromosome& a, const Chromosome& b) {
        return a.fitness > b.fitness;
    });

    // Print the best fitness and path at the end in the old population
    std::cout << "Best fitness at final iteration: " << oldPopulation[0].fitness << std::endl;
    std::cout << "Best path: ";
    for (int i = 0; i < NUM_CITIES; ++i) {
        std::cout << oldPopulation[0].path[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
