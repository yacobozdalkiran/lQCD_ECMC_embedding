//
// Created by ozdalkiran-l on 11/13/25.
//

#include <iostream>
#include "../ecmc/ecmc.h"
#include "../observables/observables.h"
#include <chrono>
#include <fstream>

void in_main_ecmc() {
    double theta_sample;
    double theta_refresh;
    int N_samples;
    double beta;
    cout << "Beta = ";
    cin >> beta;
    cout << "theta_sample = ";
    cin >> theta_sample;
    cout << "theta_refresh = ";
    cin >> theta_refresh;
    cout << "N_samples = ";
    cin >> N_samples;

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
    auto plaq = plaquette_stats(links, lat);
    cout << "Initial <P> = " << plaq.mean << " +- " << plaq.stddev << endl;

    vector<double> meas_plaquette = ecmc_samples(links, lat, beta, N_samples, theta_sample, theta_refresh, rng);

    plaq = plaquette_stats(links, lat);
    cout << "Plaquette finale : "<< endl;
    cout << "<P> = " << plaq.mean << " ± " << plaq.stddev << endl;
    string file = "plaquette_ecmc_cold.txt";
    cout << "Ecriture des données..." << endl;
    ofstream out(file);
    for (auto p : meas_plaquette) {
        out << p << " ";
    }
    out.close();

    cout << "Hot start... \n";
    hot_start(links, lat,rng);
    plaq = plaquette_stats(links, lat);
    cout << "Initial <P> = " << plaq.mean << " +- " << plaq.stddev << endl;

    meas_plaquette = ecmc_samples(links, lat, beta, N_samples, theta_sample, theta_refresh, rng);

    plaq = plaquette_stats(links, lat);
    cout << "Plaquette finale : "<< endl;
    cout << "<P> = " << plaq.mean << " ± " << plaq.stddev << endl;
    file = "plaquette_ecmc_hot.txt";
    cout << "Ecriture des données..." << endl;
    ofstream out_hot(file);
    for (auto p : meas_plaquette) {
        out_hot << p << " ";
    }
    out_hot.close();
}

int main() {
    in_main_ecmc();
    return 0;
}