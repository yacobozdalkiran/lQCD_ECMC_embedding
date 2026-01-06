cmake_minimum_required(VERSION 3.14)
project(lQCD_ECMC_embedding)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -fno-omit-frame-pointer -march=native -Wall -Wextra")

#Si Eigen est déjà installé, décommenter
find_package(Eigen3 3.3 REQUIRED NO_MODULE)


# Si Eigen n'est pas installé, décommenter les lignes suivantes : télécharge Eigen automatiquement si absent
#include(FetchContent)
#FetchContent_Declare(
#        eigen
#        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
#        GIT_TAG 3.4.0  # version stable actuelle
#)
#FetchContent_MakeAvailable(eigen)

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

target_link_libraries(test_metropolis PRIVATE Eigen3::Eigen)

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

target_link_libraries(test_ecmc PRIVATE Eigen3::Eigen)

add_executable(test_ecmc_vs_metro
        src/bins/test_ecmc_vs_metro.cpp
        src/lattice/lattice.cpp
        src/observables/observables.cpp
        src/su3utils/su3utils.cpp
        src/ecmc/ecmc.cpp
        src/metropolis/metropolis.cpp
)

target_link_libraries(test_ecmc_vs_metro PRIVATE Eigen3::Eigen)