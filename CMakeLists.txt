cmake_minimum_required(VERSION 4.0)
project(lQCD_ECMC_embedding)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -fno-omit-frame-pointer -march=native -Wall -Wextra")

add_executable(test_metropolis
        src/lattice/lattice.cpp
        src/lattice/lattice.h
        src/su3utils/su3utils.h
        src/su3utils/su3utils.cpp
        src/observables/observables.h
        src/observables/observables.cpp
        src/metropolis/metropolis.h
        src/metropolis/metropolis.cpp
        src/bins/test_metropolis.cpp)

add_executable(test_ecmc
        src/bins/test_ecmc.cpp
        src/lattice/lattice.cpp
        src/lattice/lattice.h
        src/su3utils/su3utils.h
        src/su3utils/su3utils.cpp
        src/observables/observables.h
        src/observables/observables.cpp
        src/ecmc/ecmc.h
        src/ecmc/ecmc.cpp)
