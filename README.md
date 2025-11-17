## Event-Chain and Metropolis algorithms for pure gauge SU(3) lattice QCD

This is a set of algorithms for pure gauge SU(3) lattice QCD with Wilson action. 
It contains an implementation of Event-Chain Monte Carlo with XY-Embedding and Metropolis with SU(2) embedding
to generate gauge configurations.

The compilation uses `CMake` and the package `Eigen` for linear algebra.

### Compilation

First clone the repo anywhere you want:

```bash
git clone https://github.com/yacobozdalkiran/lQCD_ECMC_embedding.git
````

Then move into the directory `lQCD_ECMC_embedding`:

 * #### If `Eigen` is installed:

```bash
mkdir build
cd build
cmake ..
make -j[NB]
```

Replace `[NB]` with the number of processors you want to use for compilation.

* #### If `Eigen` is not installed :

Go in the `CMakeLists.txt` file and uncomment this section :

```cmake
include(FetchContent)
FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0  # version stable actuelle
)
FetchContent_MakeAvailable(eigen)
```

Then comment the line 

```cmake
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
```
You can now compile :

```bash
mkdir build
cd build
cmake ..
make -j[NB]
```

### Test

If the compilation worked, you have in the `build` folder 2 executables : `test_metropolis` and `test_ecmc`.
You can run them with

```bash
./test_metropolis
```

with a result resembling 

```
Beta = 6.0
Geometry initialized : L = 4, V = 256
Hot start:  <P> = -0.00012138 ± 0.00290459

Starting Metropolis, beta = 6, epsilon = 0.15
Burn-in...
Finished ! Measuring...
Measure 0, Metropolis step 0, <P> = 0.598421 ± 0.00271545, acceptance = 0.653646
Measure 1, Metropolis step 100, <P> = 0.594656 ± 0.00277399, acceptance = 0.636719
Measure 2, Metropolis step 200, <P> = 0.589164 ± 0.00277226, acceptance = 0.630371
Measure 3, Metropolis step 300, <P> = 0.603636 ± 0.00236127, acceptance = 0.624512
Measure 4, Metropolis step 400, <P> = 0.604277 ± 0.00254303, acceptance = 0.63265
Measure 5, Metropolis step 500, <P> = 0.603004 ± 0.0026091, acceptance = 0.62207
Measure 6, Metropolis step 600, <P> = 0.60089 ± 0.00262507, acceptance = 0.63151
Measure 7, Metropolis step 700, <P> = 0.595829 ± 0.00261759, acceptance = 0.614095
Measure 8, Metropolis step 800, <P> = 0.59379 ± 0.00280752, acceptance = 0.641276
Measure 9, Metropolis step 900, <P> = 0.593844 ± 0.0025536, acceptance = 0.622884
Elapsed time : 1.14808s
Writing to file...
Done!
```

and 

```bash
./test_ecmc
```

with a result of type 

```
Beta = 6.0
theta_sample = 10000
theta_refresh = 1000
N_samples = 10
Loi poisson ? (yes : 1, no : 0) :
0
Geometry initialized : L = 4, V = 256
Cold start... 
Initial <P> = 1 +- 0
beta = 6
Sample 0, <P> = 0.726567 +- 0.00182931, 30326 events
Sample 1, <P> = 0.646757 +- 0.00215894, 29819 events
Sample 2, <P> = 0.630111 +- 0.00241089, 30376 events
Sample 3, <P> = 0.617283 +- 0.00249632, 30355 events
Sample 4, <P> = 0.615897 +- 0.00248113, 30368 events
Sample 5, <P> = 0.62518 +- 0.00227914, 30475 events
Sample 6, <P> = 0.616828 +- 0.00257435, 30374 events
Sample 7, <P> = 0.616576 +- 0.00235665, 30587 events
Sample 8, <P> = 0.611548 +- 0.00248902, 30503 events
Sample 9, <P> = 0.613121 +- 0.00246751, 30344 events
Plaquette finale : 
<P> = 0.613121 ± 0.00246751
Ecriture des données...
Hot start... 
Initial <P> = 0.00552358 +- 0.00313926
beta = 6
Sample 0, <P> = 0.371876 +- 0.00332067, 30672 events
Sample 1, <P> = 0.437061 +- 0.00333759, 32337 events
Sample 2, <P> = 0.476746 +- 0.00334412, 32461 events
Sample 3, <P> = 0.49509 +- 0.00321389, 31860 events
Sample 4, <P> = 0.518052 +- 0.00308264, 32068 events
Sample 5, <P> = 0.524819 +- 0.00307153, 32199 events
Sample 6, <P> = 0.535728 +- 0.00307301, 31709 events
Sample 7, <P> = 0.545686 +- 0.00296952, 31696 events
Sample 8, <P> = 0.560358 +- 0.00287375, 31636 events
Sample 9, <P> = 0.551432 +- 0.00300093, 31347 events
Plaquette finale : 
<P> = 0.551375 ± 0.00300084
Ecriture des données...
```