//
// Created by ozdalkiran-l on 11/17/25.
//

#include <iostream>
#include "../metropolis/metropolis.h"
#include "../ecmc/ecmc.h"
#include <chrono>
#include <fstream>

using namespace std;

void in_main_ecmc() {
    //Computes for several beta the mean plaquette with ECMC

    //Lattice params
    int L = 4;
    int Nx = L, Ny=L, Nz=L, Nt = L;
    Lattice lat(Nx, Ny, Nz, Nt);
    vector<Complex> links(lat.V*4*9);
    vector<double> beta = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0};

    //ECMC params
    double theta_sample = 1000;
    double theta_refresh = 300;
    double N_samples = 5000;
    bool poisson = false;
    double epsilon_set = 0.15;
    random_device rd;
    mt19937_64 rng(rd());

    vector<vector<double>> meas(beta.size());

    for (size_t i = 0; i < beta.size(); i++) {
        cold_start(links, lat);
        meas[i] = ecmc_samples_improved(links, lat, beta[i], N_samples, theta_sample, theta_refresh, rng, poisson, epsilon_set);
    }

    ofstream file("plaq_beta_ecmc.txt");
    for (size_t i = 0; i < beta.size(); i++) {
        file << beta[i] << " ";
        for (size_t j = 0; j < meas[i].size(); j++) {
            file << meas[i][j] << " ";
        }
        file << endl;
    }

}

void in_main_metro() {

    //Lattice
    int L = 4;
    int Nx = L, Ny=L, Nz=L, Nt = L;
    Lattice lat(Nx, Ny, Nz, Nt);
    vector<Complex> links(lat.V*4*9);

    //Params Metro
    double epsilon = 0.15;
    int n_set = 20;
    int n_meas = 100;
    int n_sweeps_meas = 20;
    int n_hits = 10;
    int n_burnin = 5000;
    size_t accepted = 0;
    size_t proposed = 0;
    vector<double> beta = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0};
    random_device rd;
    mt19937_64 rng(rd());

    vector<vector<double>> meas;
    meas.resize(beta.size());

    for (size_t i = 0; i < beta.size(); i++) {
        cold_start(links, lat);
        meas[i] = metropolis_samples(links, lat, beta[i], epsilon, n_set, n_meas, n_sweeps_meas, n_hits, n_burnin, accepted, proposed, rng);
    }

    ofstream file("plaq_beta_metro.txt");
    for (size_t i = 0; i < beta.size(); i++) {
        file << beta[i] << " ";
        for (size_t j = 0; j < meas[i].size(); j++) {
            file << meas[i][j] << " ";
        }
        file << endl;
    }
}

int main() {
    in_main_ecmc();
    return 0;
}
