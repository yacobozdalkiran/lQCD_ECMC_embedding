//
// Created by ozdalkiran-l on 11/13/25.
//

#include <iostream>
#include "../ecmc/ecmc.h"
#include "../observables/observables.h"
#include <chrono>
#include <fstream>

#include "../su3utils/su3utils.h"

void in_main_ecmc() {
    //Test du ECMC en temps continu
    int L;
    double theta_sample;
    double theta_refresh;
    int N_samples;
    double beta;
    bool poisson;
    cout << "L = ";
    cin >> L;
    cout << "Beta = ";
    cin >> beta;
    cout << "theta_sample = ";
    cin >> theta_sample;
    cout << "theta_refresh = ";
    cin >> theta_refresh;
    cout << "N_samples = ";
    cin >> N_samples;
    cout << "Poisson law ? (yes : 1, no : 0) : ";
    cin >> poisson;
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
    cout << "\nCold start... \n";
    cold_start(links, lat);
    auto plaq = plaquette_stats(links, lat);
    cout << "Initial <P> = " << plaq.mean << " +- " << plaq.stddev << endl;

    vector<double> meas_plaquette = ecmc_samples(links, lat, beta, N_samples, theta_sample, theta_refresh, rng, poisson);

    plaq = plaquette_stats(links, lat);
    cout << "Final plaquette : "<< endl;
    cout << "<P> = " << plaq.mean << " ± " << plaq.stddev << endl;
    string file = "plaquette_ecmc_cold.txt";
    cout << "Writing data..." << endl;
    ofstream out(file);
    for (auto p : meas_plaquette) {
        out << p << " ";
    }
    out.close();

    cout << "\nHot start... \n";
    hot_start(links, lat,rng);
    plaq = plaquette_stats(links, lat);
    cout << "Initial <P> = " << plaq.mean << " +- " << plaq.stddev << endl;
    cout << "Q = " << topo_charge_clover(links, lat) << endl;

    meas_plaquette = ecmc_samples(links, lat, beta, N_samples, theta_sample, theta_refresh, rng, poisson);

    plaq = plaquette_stats(links, lat);
    cout << "Final plaquette : "<< endl;
    cout << "<P> = " << plaq.mean << " ± " << plaq.stddev << endl;
    cout << "Q = " << topo_charge_clover(links, lat) << endl;
    file = "plaquette_ecmc_hot.txt";
    cout << "Writing data..." << endl;
    ofstream out_hot(file);
    for (auto p : meas_plaquette) {
        out_hot << p << " ";
    }
    out_hot.close();
}

void in_main_ecmc_args(int L, double beta, double theta_sample, double theta_refresh, int N_samples, int poisson) {
    //Test du ECMC en temps continu
    cout << "L = " << L << std::endl;
    cout << "Beta = " << beta << std::endl;
    cout << "theta_sample = " << theta_sample << std::endl;
    cout << "theta_refresh = " << theta_refresh << std::endl;
    cout << "N_samples = " << N_samples << std::endl;
    cout << "Poisson law ? (yes : 1, no : 0) : " << poisson << std::endl;
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
    cout << "\nCold start... \n";
    cold_start(links, lat);
    auto plaq = plaquette_stats(links, lat);
    cout << "Initial <P> = " << plaq.mean << " +- " << plaq.stddev << endl;

    vector<double> meas_plaquette = ecmc_samples(links, lat, beta, N_samples, theta_sample, theta_refresh, rng, poisson);

    plaq = plaquette_stats(links, lat);
    cout << "Final plaquette : "<< endl;
    cout << "<P> = " << plaq.mean << " ± " << plaq.stddev << endl;
    string file = "plaquette_ecmc_cold.txt";
    cout << "Writing data..." << endl;
    ofstream out(file);
    for (auto p : meas_plaquette) {
        out << p << " ";
    }
    out.close();

    cout << "\nHot start... \n";
    hot_start(links, lat,rng);
    plaq = plaquette_stats(links, lat);
    cout << "Initial <P> = " << plaq.mean << " +- " << plaq.stddev << endl;
    cout << "Q = " << topo_charge_clover(links, lat) << endl;

    meas_plaquette = ecmc_samples(links, lat, beta, N_samples, theta_sample, theta_refresh, rng, poisson);

    plaq = plaquette_stats(links, lat);
    cout << "Final plaquette : "<< endl;
    cout << "<P> = " << plaq.mean << " ± " << plaq.stddev << endl;
    cout << "Q = " << topo_charge_clover(links, lat) << endl;
    file = "plaquette_ecmc_hot.txt";
    cout << "Writing data..." << endl;
    ofstream out_hot(file);
    for (auto p : meas_plaquette) {
        out_hot << p << " ";
    }
    out_hot.close();
}


void in_main_ecmc_improved(double epsilon_set) {
    //Test du ECMC en temps continu
    int L;
    double theta_sample;
    double theta_refresh;
    int N_samples;
    double beta;
    bool poisson;
    cout << "L = ";
    cin >> L;
    cout << "Beta = ";
    cin >> beta;
    cout << "theta_sample = ";
    cin >> theta_sample;
    cout << "theta_refresh = ";
    cin >> theta_refresh;
    cout << "N_samples = ";
    cin >> N_samples;
    cout << "Poisson law ? (yes : 1, no : 0) : ";
    cin >> poisson;
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
    //cout << "\nCold start... \n";
    hot_start(links, lat, rng);
    auto plaq = plaquette_stats(links, lat);
    cout << "Initial <P> = " << plaq.mean << " +- " << plaq.stddev << endl;
    cout << "Inital S = " << wilson_action(links, lat, beta) << endl;

    auto start = chrono::high_resolution_clock::now();
    vector<double> meas_plaquette = ecmc_samples_improved(links, lat, beta, N_samples, theta_sample, theta_refresh, rng, poisson,epsilon_set);
    auto end = chrono::high_resolution_clock::now();
    auto elapsed = end - start;
    cout << "Elapsed time : " << chrono::duration_cast<chrono::seconds>(elapsed).count() << " s" << endl;

    plaq = plaquette_stats(links, lat);
    cout << "Final plaquette : "<< endl;
    cout << "<P> = " << plaq.mean << " ± " << plaq.stddev << endl;
    cout << "Final S = " << wilson_action(links, lat, beta) << endl;
    string file = "plaquette_ecmci_cold.txt";
    cout << "Writing data..." << endl;
    ofstream out(file);
    for (auto p : meas_plaquette) {
        out << p << " ";
    }
    out.close();

    cout << "\nHot start... \n";
    hot_start(links, lat,rng);
    plaq = plaquette_stats(links, lat);
    cout << "Initial <P> = " << plaq.mean << " +- " << plaq.stddev << endl;
    cout << "Q = " << topo_charge_clover(links, lat) << endl;

    meas_plaquette = ecmc_samples_improved(links, lat, beta, N_samples, theta_sample, theta_refresh, rng, poisson, epsilon_set);

    plaq = plaquette_stats(links, lat);
    cout << "Final plaquette : "<< endl;
    cout << "<P> = " << plaq.mean << " ± " << plaq.stddev << endl;
    cout << "Q = " << topo_charge_clover(links, lat) << endl;
    file = "plaquette_ecmci_hot.txt";
    cout << "Writing data..." << endl;
    ofstream out_hot(file);
    for (auto p : meas_plaquette) {
        out_hot << p << " ";
    }
    out_hot.close();
}


void in_main_ecmc_improved_args(int L, double beta, double theta_sample, double theta_refresh, int N_samples, int poisson, double epsilon_set) {
    //Test du ECMC en temps continu
    cout << "L = " << L << std::endl;
    cout << "Beta = " << beta << std::endl;
    cout << "theta_sample = " << theta_sample << std::endl;
    cout << "theta_refresh = " << theta_refresh << std::endl;
    cout << "N_samples = " << N_samples << std::endl;
    cout << "Poisson law ? (yes : 1, no : 0) : " << poisson << std::endl;
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
    //cout << "\nCold start... \n";
    hot_start(links, lat, rng);
    auto plaq = plaquette_stats(links, lat);
    cout << "Initial <P> = " << plaq.mean << " +- " << plaq.stddev << endl;
    cout << "Inital S = " << wilson_action(links, lat, beta) << endl;

    auto start = chrono::high_resolution_clock::now();
    vector<double> meas_plaquette = ecmc_samples_improved(links, lat, beta, N_samples, theta_sample, theta_refresh, rng, poisson,epsilon_set);
    auto end = chrono::high_resolution_clock::now();
    auto elapsed = end - start;
    cout << "Elapsed time : " << chrono::duration_cast<chrono::seconds>(elapsed).count() << " s" << endl;

    plaq = plaquette_stats(links, lat);
    cout << "Final plaquette : "<< endl;
    cout << "<P> = " << plaq.mean << " ± " << plaq.stddev << endl;
    cout << "Final S = " << wilson_action(links, lat, beta) << endl;
    string file = "plaquette_ecmci_cold.txt";
    cout << "Writing data..." << endl;
    ofstream out(file);
    for (auto p : meas_plaquette) {
        out << p << " ";
    }
    out.close();

    cout << "\nHot start... \n";
    hot_start(links, lat,rng);
    plaq = plaquette_stats(links, lat);
    cout << "Initial <P> = " << plaq.mean << " +- " << plaq.stddev << endl;
    cout << "Q = " << topo_charge_clover(links, lat) << endl;

    meas_plaquette = ecmc_samples_improved(links, lat, beta, N_samples, theta_sample, theta_refresh, rng, poisson, epsilon_set);

    plaq = plaquette_stats(links, lat);
    cout << "Final plaquette : "<< endl;
    cout << "<P> = " << plaq.mean << " ± " << plaq.stddev << endl;
    cout << "Q = " << topo_charge_clover(links, lat) << endl;
    file = "plaquette_ecmci_hot.txt";
    cout << "Writing data..." << endl;
    ofstream out_hot(file);
    for (auto p : meas_plaquette) {
        out_hot << p << " ";
    }
    out_hot.close();
}


int main(int argc, char* argv[]) {
    // bool improved;
    // cout << "Improved ? (1:yes/0:no) : ";
    // cin >> improved;
    // if (improved) {
    //     double epsilon_set;
    //     cout << "epsilon_set = ";
    //     cin >> epsilon_set;
    //     in_main_ecmc_improved(epsilon_set);
    // }
    // else in_main_ecmc();
    //in_main_ecmc_improved(0.15);
    if (argc != 7) cerr << "Wrong number of args, usage : <L> <beta> <theta_sample> <theta_refresh> <N_samples> <poisson>" << std::endl;
    in_main_ecmc_args(atoi(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]), atoi(argv[5]), atoi(argv[6]));
    return 0;
}