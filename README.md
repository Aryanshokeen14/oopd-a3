# CSE600 Assignment 3 – Parallel Pandemic Simulation

**Author:** Aryan Shokeen  
**Student ID:** MT26116

## Overview
Extends the sequential simulation with multithreaded execution.

## Implemented Features
- Part A sequential baseline.
- Thread-per-day implementation.
- Thread pool implementation.
- Row-strip partitioning.
- Square-block partitioning.
- Interleaved-row partitioning.
- Correctness comparison against sequential execution.
- Testing with multiple thread counts.

## Technologies
- C++17
- std::thread
- Mutexes and condition variables
- Futures and thread pool

## Compilation and Execution
```bash
make
./pandemic_sim
```

## Author
Aryan Shokeen – MT26116