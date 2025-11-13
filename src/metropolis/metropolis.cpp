//
// Created by ozdalkiran-l on 11/13/25.
//

#include "metropolis.h"
#include <chrono>
#include "../observables/observables.h"

void metropolis_sweep(vector<Complex> &links, const Lattice &lat, double beta, mt19937_64 &rng, vector<SU3> &set, size_t &accepted, size_t &proposed, int n_hits) {
    //Effectue un sweep metropolis avec n_hits hits à chaque lien.
    SU3 staple;
    SU3 Uold;
    SU3 Unew;
    accepted = 0;
    proposed = 0;
    int i_set = 0;
    uniform_int_distribution<int> index_set(0, set.size() - 1);
    uniform_real_distribution<double> unif(0.0, 1.0);
    for (size_t site = 0; site < lat.V; site++) {
        for (int mu = 0; mu < 4; mu++) {

            //On copie la staple et Uold
            compute_staple(links, lat, site, mu, staple);
            auto Umap = view_link(links, site, mu);
            Uold = Umap;

            for (int i =0; i < n_hits; i++) {
                //On choisit une matrice d'update dans set
                i_set = index_set(rng);

                //On copie Unew
                Unew = set[i_set] * Uold;

                //On propose le move
                double old_tr = (Uold * staple).trace().real();
                double new_tr = (Unew * staple).trace().real();
                double dS = -(beta/3.0) * (new_tr - old_tr);

                ++proposed;
                bool accept = false;
                if ((dS <= 0.0)||(unif(rng) < exp(-dS))) accept = true;

                if (accept) {
                    Umap = Unew;
                    ++accepted;
                }
            }
        }
    }
}

vector<double> metropolis_samples(vector<Complex> &links, const Lattice &lat, double beta, double epsilon, int n_set, int n_meas, int n_sweeps_meas, int n_hits, int n_burnin, size_t &accepted, size_t &proposed, mt19937_64 &rng) {
    //Effectue l'algo metropolis
    cout << "Création du set de matrices...\n";
    int size_set = 50;
    auto set = metropolis_set(epsilon,size_set,rng);

    cout << "\nDémarrage du Metropolis, beta = " << beta <<"\n";
    cout << "Burn-in...\n";
    for (int i = 0; i < n_burnin; i++) {
        if (i%n_set ==0) set = metropolis_set(epsilon,size_set,rng);
        metropolis_sweep(links, lat, beta, rng, set, accepted, proposed, n_hits);
    }
    cout << "Finished ! Measuring...\n";
    vector<double> measures(n_meas);
    int n_sweeps = n_meas*n_sweeps_meas;
    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < n_sweeps; i++) {
        if (i%n_set ==0) set = metropolis_set(epsilon,size_set,rng);
        metropolis_sweep(links, lat, beta, rng, set, accepted, proposed, n_hits);
        if (i%n_sweeps_meas==0) {
            auto plaq = plaquette_stats(links, lat);
            measures[i/n_sweeps_meas] = plaq.mean;
            cout << "Measure " << i/n_sweeps_meas << ", " << "Metropolis step "<< i <<", <P> = " << plaq.mean << " ± " << plaq.stddev << ", acceptance = " << double(accepted) / double(proposed) << endl;
        }
    }
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Elapsed time : " << elapsed.count() << "s" << endl;
    return measures;
}