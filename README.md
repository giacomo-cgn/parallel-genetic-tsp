# Parallel Genetic TSP project
This project is a parallelized genetic algorithm for solving the Travelling Salesman Problem in C++ done for the SPM course of the University of Pisa.

The genetic algorithm employs crossover and random mutation to evolve. The proposed path to solve TSP is considered as theset of genes of each individual. 
Two different parallel implementations of the algorithm are provided, togehter with the sequential version
Tested on [Zimbabwe - 929 Cities](https://www.math.uwaterloo.ca/tsp/world/countries.html) dataset of cities.

`Report_SPM.pdf` contains more detail on the implementation and experimental results of the project.

`plots.ipynb` is a notebook for plotting useful graphs on experiments results.

The code is contained in the `/src` directory.

## Genetic algorithm
It consists of the following steps:
1. Generation of the adjacency matrix of city distances 
2. Random initialization of parent and child populations
3. Iterate for _n_ steps:
   - Sort individuals by fitness (TSP path length)
   - Select the top _elitism_rate_ individuals for crossover
   - Crossover to generate children (same amount of individuals as parents)
   - Mutate for each child between 0 and _num_cities_ $\times$ _mutation_ratio_ genes
   - Calcuate fitness for children
  
## Parallelization
It is a data parallel context, as each iteration of the algorithm cannot begin until the previous has ended.
Each phase of the algorith was parallelized with a _map_ parallel pattern, following a fork-join model.

The first implementation only uses C++ native threads. The `ParallelMap` class implements the _map_ pattern.

In the second implementation, the construct `parallel_for` from [FastFlow](https://github.com/fastflow/fastflow) library was used instead.

## Experiments
Experiments on performance were run on a 32 core machine. The parallel versions of the program have been tested with varying number of workers and confronted with the sequential version.
Speedup is the main measure of performance.

We attained a maximum speedup of 5.30 with 10 workers for the thread implementation and 6.98 with 14 workers for the FastFlow implementation.

## How to execute
Open Terminal:
Open the terminal on your Ubuntu system.

- Inside `/parallel-genetic-tsp` create a build directory:
```
mkdir build
cd build
```
- Generate Makefiles with CMake:
```
cmake ..
```
The `..` points to the parent directory where CMakeLists.txt is located. Adjust the path accordingly if CMakeLists.txt is in a different location.

Build the program:
```
make
```

Run the program:
```
./ParallelGeneticTSP
```

All experiment will be runned in wxecuting in this way.
If you want to explore different settings for experiments or only reproduce part of it, modify `/src/main.cpp` and run build the program again with `cmake`.



   

