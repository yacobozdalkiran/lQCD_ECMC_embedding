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
    int L = 4;
    int Nx = L, Ny = L, Nz = L, Nt = L;
    random_device rd;
    mt19937_64 rng(rd());
    //On construit la lattice de sites
    Lattice lat(Nx, Ny, Nz, Nt);
    size_t V = lat.V;

    cout << "Geometry initialized : L = " << L << ", V = " << V << "\n";

    //On construit le vector de liens en mémoire contigue
    vector<Complex> links(4*V*9);

    //On initialise tout à l'identité
    cout << "Cold start... \n";
    cold_start(links, lat);
    auto cold = plaquette_stats(links, lat);
    cout << "Cold start:  <P> = " << cold.mean << " ± " << cold.stddev << endl;

    //Hot start
    cout << "Hot start... \n";
    hot_start(links, lat,rng);
    auto hot = plaquette_stats(links, lat);
    cout << "Hot start:  <P> = " << hot.mean << " ± " << hot.stddev << endl;

    cout << "Création du set de matrices...\n";
    double epsilon = 0.15;
    int size_set = 50;
    cout << "epsilon = " << epsilon << "\n";
    auto set = metropolis_set(epsilon,size_set,rng);
    //cout << "Déterminant de la première : " << set[1].determinant() << endl;
    //cout << "Distance à l'identité de la première : " << sqrt(((set[2]-SU3::Identity())*(set[1] - SU3::Identity())).trace().real()) << endl;


    cout << "\nDémarrage du Metropolis, beta = " << beta <<"\n";
    size_t accepted = 0;
    size_t proposed = 0;
    int n_hits = 6;
    int n_burnin = 2000; //Burn-in a 2000 pour L=4 beta=6
    int n_sweeps = 5000;
    int n_set = 20; //Refresh du set tous les n_set sweeps
    int n_meas = 100; //Mesure tous les n_meas sweeps

    cout << "Burn-in...\n";
    for (int i = 0; i < n_burnin; i++) {
        if (i%n_set ==0) set = metropolis_set(epsilon,size_set,rng);
        metropolis_sweep(links, lat, beta, rng, set, accepted, proposed, n_hits);
    }
    cout << "Finished ! Measuring...\n";

    vector<double> measures(n_sweeps/n_meas);

    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < n_sweeps; i++) {
        if (i%n_set ==0) set = metropolis_set(epsilon,size_set,rng);
        metropolis_sweep(links, lat, beta, rng, set, accepted, proposed, n_hits);
        if (i%n_meas==0) {
            auto plaq = plaquette_stats(links, lat);
            measures[i/n_meas] = plaq.mean;
            cout << "Measure " << i/n_meas << ", " << "Metropolis step "<< i <<", <P> = " << plaq.mean << " ± " << plaq.stddev << ", acceptance = " << double(accepted) / double(proposed) << endl;
        }
    }
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Elapsed time : " << elapsed.count() << "s" << endl;
    cout << "Writing to file...\n";
    ofstream file("metro_plaquette_hot.txt");
    if (!file) cerr << "Can't open file" << endl;
    for (int i = 0; i < n_sweeps/n_meas; i++) {
        file << measures[i] << " ";
    }
    file.close();
    cout << "Done!\n";
}

int main() {
    in_main_metropolis();
    return 0;
}