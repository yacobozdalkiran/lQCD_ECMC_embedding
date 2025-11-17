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
Poisson law ? (yes : 1, no : 0) : 0
Geometry initialized : L = 4, V = 256

Cold start... 
Initial <P> = 1 +- 0
beta = 6
Sample 0, <P> = 0.736834 +- 0.00178981, 30491 events
Sample 1, <P> = 0.657786 +- 0.00225725, 29722 events
Sample 2, <P> = 0.627181 +- 0.00240389, 30013 events
Sample 3, <P> = 0.62207 +- 0.00246158, 30260 events
Sample 4, <P> = 0.618506 +- 0.00230837, 30316 events
Sample 5, <P> = 0.615442 +- 0.00245077, 30306 events
Sample 6, <P> = 0.604019 +- 0.00240745, 30507 events
Sample 7, <P> = 0.613388 +- 0.00248805, 30474 events
Sample 8, <P> = 0.613853 +- 0.0023788, 30361 events
Sample 9, <P> = 0.612192 +- 0.00257995, 30518 events
Final plaquette : 
<P> = 0.612192 ± 0.00257995
Writing data...

Hot start... 
Initial <P> = -0.00809392 +- 0.00287814
beta = 6
Sample 0, <P> = 0.382902 +- 0.00331101, 30963 events
Sample 1, <P> = 0.442582 +- 0.00327299, 32615 events
Sample 2, <P> = 0.48202 +- 0.00314621, 32843 events
Sample 3, <P> = 0.497703 +- 0.00327575, 32131 events
Sample 4, <P> = 0.513254 +- 0.00311604, 31757 events
Sample 5, <P> = 0.526483 +- 0.00299539, 31708 events
Sample 6, <P> = 0.53425 +- 0.00290803, 31860 events
Sample 7, <P> = 0.539642 +- 0.00295148, 31675 events
Sample 8, <P> = 0.541349 +- 0.00292008, 31739 events
Sample 9, <P> = 0.541747 +- 0.00306963, 31548 events
Final plaquette : 
<P> = 0.541736 ± 0.00306955
Writing data...

```