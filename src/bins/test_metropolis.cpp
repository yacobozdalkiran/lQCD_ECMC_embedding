//
// Created by ozdalkiran-l on 11/13/25.
//

#include <iostream>
#include "../metropolis/metropolis.h"
#include "../observables/observables.h"
#include <chrono>
#include <fstream>

void in_main_metropolis() {
    double beta;
    cout << "Beta = ";
    cin >> beta;
    int L = 6;
    int Nx = L, Ny = L, Nz = L, Nt = L;
    random_device rd;
    mt19937_64 rng(rd());
    //On construit la lattice de sites
    Lattice lat(Nx, Ny, Nz, Nt);
    size_t V = lat.V;

    cout << "Geometry initialized : L = " << L << ", V = " << V << "\n";

    //On construit le vector de liens de jauge en mémoire contigue
    vector<Complex> links(4*V*9);

    //Hot start
    hot_start(links, lat,rng);
    auto hot = plaquette_stats(links, lat);
    cout << "Hot start:  <P> = " << hot.mean << " ± " << hot.stddev << endl;

    double epsilon = 0.15;
    int n_set = 20; //Refresh du set tous les n_set sweeps
    int n_meas = 10; //n_meas mesures de plaquettes
    int n_sweeps_meas = 100; //n_sweeps_meas sweeps entre chaque mesure
    int n_hits = 6; // n_hits hits par lien pour chaque sweep
    int n_burnin = 2000; //Burn-in a 2000 pour L=4 beta=6

    size_t accepted = 0;
    size_t proposed = 0;

    vector<double> measures = metropolis_samples(links, lat, beta, epsilon, n_set, n_meas, n_sweeps_meas, n_hits, n_burnin, accepted, proposed, rng);

    cout << "Writing to file...\n";
    ofstream file("metro_plaquette_hot.txt");
    if (!file) cerr << "Can't open file" << endl;
    for (const auto &mea : measures) {
        file << mea << " ";
    }
    file.close();
    cout << "Done!\n";
}

int main() {
    in_main_metropolis();
    return 0;
}