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
    double theta_sample;
    double theta_refresh;
    int N_samples;
    double beta;
    bool poisson;
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

    meas_plaquette = ecmc_samples(links, lat, beta, N_samples, theta_sample, theta_refresh, rng, poisson);

    plaq = plaquette_stats(links, lat);
    cout << "Final plaquette : "<< endl;
    cout << "<P> = " << plaq.mean << " ± " << plaq.stddev << endl;
    file = "plaquette_ecmc_hot.txt";
    cout << "Writing data..." << endl;
    ofstream out_hot(file);
    for (auto p : meas_plaquette) {
        out_hot << p << " ";
    }
    out_hot.close();
}

void in_main_ecmc_d() {
    //Test du ECMC en temps discret
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

    //ECMC temps discret, initialisation
    uniform_int_distribution<size_t> site_dist(0, lat.V-1);
    uniform_int_distribution<int> mu_dist(0,3);
    size_t site = site_dist(rng);
    int epsilon = 1;
    int mu = mu_dist(rng);
    array<SU3, 6> list_staple;
    SU3 R = random_su3(rng);
    double beta = 6.0;
    double eta = 0.5;
    SU3 lambda_3;
    lambda_3 << Complex(1.0,0.0), Complex(0.0,0.0), Complex(0.0,0.0),
                Complex(0.0,0.0), Complex(-1.0,0.0), Complex(0.0, 0.0),
                Complex(0.0,0.0), Complex(0.0,0.0), Complex(0.0, 0.0);

    for (int i = 0; i<100000; i++) {
        compute_list_staples(links, lat, site, mu, list_staple); //On calcule les 6 staples
        int j = update_until_reject_d(links, site, mu, list_staple, R, epsilon, beta, eta, rng); //On update le lien jusqu'à un rejet par la plaquette j
        auto l = lift(links, lat, site, mu, j, R, lambda_3, rng); //On lifte dans la plaquette j
        site = l.first.first;
        mu = l.first.second;
        epsilon = l.second;
        if (i%1000 == 0) {
            //Les mesures sont faites lors d'un event -> biais, mais semble similaire à l'algo en temps continu
            plaq = plaquette_stats(links, lat);
            cout << "<P> = " << plaq.mean << " +- " << plaq.stddev << endl;
            //On refresh la chaîne (le lien et R)
            site = site_dist(rng);
            epsilon = 1;
            mu = mu_dist(rng);
            R = random_su3(rng);
        }
    }


}

int main() {
    in_main_ecmc();
    return 0;
}