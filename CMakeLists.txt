cmake_minimum_required(VERSION 3.16)
project(ParallelGeneticTSP)

include(FetchContent)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE -O3)

find_package(Threads REQUIRED)

# Include FastFlow library
message(STATUS "Fetching 'FastFlow' from GitHub...")
FetchContent_Declare(
        fastflow
        GIT_REPOSITORY https://github.com/fastflow/fastflow.git
)
FetchContent_MakeAvailable(fastflow)

# !!! Delete if causing error in CMakelists.txt of FastFlow if "pthread" error arises

# FastFlow is header-only
include_directories(${fastflow_SOURCE_DIR})

# include_directories(include)



# Add source files for the executable
add_executable(ParallelGeneticTSP
    src/utimer.h
    src/utils.h
    src/main.cpp 
    src/sequential.cpp
    src/threads.cpp
    src/utils.cpp
    )

# Link FastFlow library to the executable (if needed)
# target_link_libraries(ParallelGeneticTSP PRIVATE FastFlow::FastFlow)
