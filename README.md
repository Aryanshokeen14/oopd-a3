# CSE600 Assignment 3 – Asynchronous Pandemic Simulation

**Author:** Aryan Shokeen  
**Student ID:** MT26116

## Overview
Extends the parallel simulation with asynchronous snapshot writing and parameter sweeps.

## Implemented Features
- All Part A and Part B features.
- Blocking PPM snapshot generation.
- Asynchronous PPM snapshot generation.
- Configurable snapshot interval.
- Maximum pending asynchronous writes.
- Asynchronous parameter sweeps.
- Worker limits.
- Time-limited execution.
- Completion-order result collection.

## Technologies
- C++17
- std::async
- std::future
- Multithreading

## Compilation and Execution
```bash
make
./pandemic_sim
```

## Author
Aryan Shokeen – MT26116