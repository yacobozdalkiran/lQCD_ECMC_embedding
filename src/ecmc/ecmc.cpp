//
// Created by ozdalkiran-l on 11/13/25.
//

#include "ecmc.h"
#include "../lattice/lattice.h"
#include "../su3utils/su3utils.h"
#include "../observables/observables.h"

void compute_list_staples(const vector<Complex> &links, const Lattice &lat, size_t site, int mu, array<SU3,6> &list_staple) {
    //Calcule la liste des staples autour un lien de jauge
    int index = 0;
    for (int nu = 0; nu < 4; nu++) {
        if (nu == mu) {
            continue;
        }
        size_t x = site; //x
        size_t xmu = lat.neighbor[x][mu][0]; //x+mu
        size_t xnu = lat.neighbor[x][nu][0]; //x+nu
        size_t xmunu = lat.neighbor[xmu][nu][1]; //x+mu-nu
        size_t xmnu = lat.neighbor[x][nu][1]; //x-nu
        auto U0 = view_link_const(links, xmu, nu);
        auto U1 = view_link_const(links, xnu, mu);
        auto U2 = view_link_const(links, x, nu);
        list_staple[index] = U0 * U1.adjoint() * U2.adjoint();
        auto V0 = view_link_const(links, xmunu, nu);
        auto V1 = view_link_const(links, xmnu, mu);
        auto V2 = view_link_const(links, xmnu, nu);
        list_staple[index+1] = V0.adjoint() * V1.adjoint() * V2;
        index += 2;
    }
}

void compute_reject(double A, double B, double &gamma, double &reject, int epsilon) {
    //TODO: A refactoriser pour rendre plus clair
    if (epsilon == -1) B = -B;
    double R = sqrt(A*A + B*B);
    double phi = atan2(-A/R, B/R);

    if (phi<0) phi += 2*M_PI;
    //cout << "phi = " << phi << endl;
    double period = 0.0, p1 = 0.0, p2 = 0.0;
    array<double,4> intervals = {0.0, 0.0, 2*M_PI, 2*M_PI};
    if (phi < M_PI/2.0) {
        //cout << "cas 1"<< endl;
        intervals[1] = M_PI/2.0 + phi;
        intervals[2] = 3*M_PI/2.0 + phi;
        p1 = R * ( sin(intervals[1] - phi) - sin(intervals[0] - phi) );
        p2 = R * ( sin(intervals[3] - phi) - sin(intervals[2] - phi) );
        if ((p1<0)&&(p2<0)) cerr << "Périodes négatives !"<< endl;
        period = p1+p2;
        //cout << "contrib periodique = " << period << endl;
        gamma = gamma - std::floor(gamma/period)*period;
        if (gamma>p1) {
            gamma -= p1;
            double alpha = gamma/R + sin(intervals[2]-phi);
            double theta1 = fmod((phi + asin(alpha) + 2* M_PI), 2*M_PI);
            double theta2 = fmod((phi + M_PI - asin(alpha) + 2* M_PI), 2*M_PI);
            if ((theta1<intervals[3])&&(theta1>intervals[2])){
                reject = theta1;
            }
            else {
                reject = theta2;
            }
        }
        else {
            double alpha = gamma/R + sin(intervals[0]-phi);
            double theta1 = fmod((phi + asin(alpha) + 2* M_PI), 2*M_PI);
            double theta2 = fmod((phi + M_PI - asin(alpha) + 2* M_PI), 2*M_PI);
            if ((theta1<intervals[1])&&(theta1>intervals[0])){
                reject = theta1;
            }
            else {
                reject = theta2;
            }
        }
    }
    if (phi > 3*M_PI/2.0) {
        //cout << "cas 2" << endl;
        intervals[1] = -3*M_PI/2.0 + phi;
        intervals[2] = -M_PI/2.0 + phi;
        //cout << "[" << intervals[0] << ", " << intervals[1] << "]" << endl;
        //cout << "[" << intervals[2] << ", " << intervals[3] << "]" << endl;
        p1 = R * ( sin(intervals[1] - phi) - sin(intervals[0] - phi) );
        p2 = R * ( sin(intervals[3] - phi) - sin(intervals[2] - phi) );
        if ((p1<0)&&(p2<0)) cerr << "Périodes négatives !"<< endl;
        period = p1+p2;
        //cout << "contrib periodique = " << period << endl;
        gamma = gamma - std::floor(gamma/period)*period;
        if (gamma>p1) {
            gamma -= p1;
            double alpha = gamma/R + sin(intervals[2]-phi);
            double theta1 = fmod((phi + asin(alpha) + 2* M_PI), 2*M_PI);
            double theta2 = fmod((phi + M_PI - asin(alpha) + 2* M_PI), 2*M_PI);
            if ((theta1<intervals[3])&&(theta1>intervals[2])){
                reject = theta1;
            }
            else {
                reject = theta2;
            }
        }
        else {
            double alpha = gamma/R + sin(intervals[0]-phi);
            double theta1 = fmod((phi + asin(alpha) + 2* M_PI), 2*M_PI);
            double theta2 = fmod((phi + M_PI - asin(alpha) + 2* M_PI), 2*M_PI);
            if ((theta1<intervals[1])&&(theta1>intervals[0])){
                reject = theta1;
            }
            else {
                reject = theta2;
            }
        }
    }
    if ((phi >= M_PI/2.0)&&(phi<= 3*M_PI/2.0)) {
        //cout << "cas 3" << endl;
        intervals[0] = -M_PI/2.0 + phi;
        intervals[1] = M_PI/2.0 + phi;
        period = R * ( sin(intervals[1] - phi) - sin(intervals[0] - phi) );
        if (period<0) cerr << "Période négative !"<< endl;
        //cout << "contrib periodique = " << period << endl;
        gamma = gamma - std::floor(gamma/period)*period;
        double alpha = gamma/R + sin(intervals[0]-phi);
        double theta1 = fmod((phi + asin(alpha) + 2* M_PI), 2*M_PI);
        double theta2 = fmod((phi + M_PI - asin(alpha) + 2* M_PI), 2*M_PI);
        if ((theta1<intervals[1])&&(theta1>intervals[0])){
            reject = theta1;
        }
        else {
            reject = theta2;
        }
    }
}

void compute_reject_angles(const vector<Complex> &links, size_t site, int mu, const array<SU3,6> &list_staple, const SU3 &R, int epsilon, const double &beta, array<double,6> &reject_angles, mt19937_64 &rng) {
    //Calcule la liste des 6 angles de rejets pour le lien (site, mu) (1 angle par plaquette associée au lien)
    double gamma;
    SU3 P;
    double A = 0;
    double B = 0;

    uniform_real_distribution<double> unif(0.0,1.0);
    for (int i = 0; i < 6; i++) {
        gamma = -log(unif(rng));

        P = R.adjoint() * view_link_const(links, site, mu) * list_staple[i] * R;
        A = P(0,0).real() + P(1,1).real();
        B = -P(0,0).imag() + P(1,1).imag();
        A *= -(beta/3.0);
        B *= -(beta/3.0);

        compute_reject(A, B, gamma, reject_angles[i], epsilon);
    }
}

int selectVariable(const vector<double> &probas, mt19937_64 &rng) {
    //Choisit un index entre 0 et probas.size()-1 selon la méthode tower of probas
    uniform_real_distribution<double> unif(0.0,1.0);
    double r = unif(rng);
    double s = 0.0;
    for (int i = 0; i < probas.size(); i++) {
        s += probas[i];
        if (s>r) {
           return i;
        }
    }
    cerr << "SelectVariable Error" << endl;
    return -1;
}


pair<pair<size_t, int>,int> lift(const vector<Complex> &links, const Lattice &lat, size_t site, int mu, int j, const SU3 &R, const SU3 &lambda_3, mt19937_64 &rng) {
    SU3 U0 = view_link_const(links, site, mu);
    auto links_staple_j = lat.staples[site][mu][j]; //Liste des 3 liens de la staple j
    SU3 U1 = view_link_const(links, links_staple_j[0].first, links_staple_j[0].second); //Les 3 matrices SU3 associées
    SU3 U2 = view_link_const(links, links_staple_j[1].first, links_staple_j[1].second);
    SU3 U3 = view_link_const(links, links_staple_j[2].first, links_staple_j[2].second);

    array<pair<size_t, int>,4> links_plaquette_j; //On rajoute le lien actuel
    links_plaquette_j[0] = make_pair(site, mu);
    links_plaquette_j[1] = links_staple_j[0];
    links_plaquette_j[2] = links_staple_j[1];
    links_plaquette_j[3] = links_staple_j[2];
    //cout << "List des sites de la plaquette j :" << endl;
    //cout << links_plaquette_j[0].first << endl << links_plaquette_j[1].first << endl << links_plaquette_j[2].first << endl << links_plaquette_j[3].first << endl;
    //cout << "List des mu de la plaquette j :" << endl;
    //cout << links_plaquette_j[0].second << endl << links_plaquette_j[1].second << endl << links_plaquette_j[2].second << endl << links_plaquette_j[3].second << endl;
    vector<double> probas(4);
    double sum =0.0;
    vector<int> sign_dS(4);
    vector<SU3> P(4);

    if (j%2 == 0) { //Forward plaquette
        P[0] = U0 * U1 * U2.adjoint() * U3.adjoint();
        P[1] = U1 * U2.adjoint() * U3.adjoint() * U0;
        P[2] = U2 * U1.adjoint() * U0.adjoint()* U3;
        P[3] = U3 * U2 * U1.adjoint() * U0.adjoint();
    }
    else { //Backward plaquette
        P[0] = U0 * U1.adjoint() * U2.adjoint() * U3;
        P[1] = U1 * U0.adjoint() * U3.adjoint() * U2;
        P[2] = U2 * U1 * U0.adjoint()* U3.adjoint();
        P[3] = U3 * U0 * U1.adjoint() * U2.adjoint();
    }
    for (int i = 0; i < 4; i++) {
        probas[i] = -(Complex(0.0,1.0) * lambda_3 * R.adjoint() * P[i]*R).trace().real();
        sign_dS[i] = dsign(probas[i]);
        probas[i] = abs(probas[i]);
        sum += probas[i];
    }
    for (int i = 0; i < 4; i++) {
        probas[i] /= sum;
    }
    //cout << "probas = [" << probas[0] << ", " << probas[1] << ", " << probas[2] << ", " << probas[3] << "]" << endl;
    //cout << "sum probas = " << probas[0] + probas[1] + probas[2] + probas[3] << endl;
    //cout << "sum = " << sum << endl;
    //cout << "signs = [" << sign_dS[0] << ", " << sign_dS[1] << ", " << sign_dS[2] << ", " << sign_dS[3] << "]" << endl;
    int index_lift = selectVariable(probas, rng);
    //cout << "index_lift = " << index_lift << endl;
    return make_pair(links_plaquette_j[index_lift], -sign_dS[index_lift]);
}

void ecmc_update(vector<Complex> &links, size_t site, int mu, double theta, int epsilon, const SU3 &R) {
    SU3 Uold = view_link_const(links, site, mu);
    view_link(links, site, mu) = R*el_3(epsilon*theta)*R.adjoint()*Uold;
    //projection_su3(links, site, mu);
}

vector<double> ecmc_samples(vector<Complex> &links, const Lattice &lat, double beta, int N_samples, double param_theta_sample, double param_theta_refresh, mt19937_64 &rng, bool poisson) {
    size_t V = lat.V;

    //Variables aléatoires
    uniform_int_distribution<size_t> random_site(0, V-1);
    uniform_int_distribution<int> random_dir(0,3);
    uniform_int_distribution<int> random_eps(0,1);
    exponential_distribution<double> random_theta_sample(1.0/param_theta_sample);
    exponential_distribution<double> random_theta_refresh(1.0/param_theta_refresh);


    //Matrice lambda_3 de Gell-Mann
    SU3 lambda_3;
    lambda_3 << Complex(1.0,0.0), Complex(0.0,0.0), Complex(0.0,0.0),
                Complex(0.0,0.0), Complex(-1.0,0.0), Complex(0.0, 0.0),
                Complex(0.0,0.0), Complex(0.0,0.0), Complex(0.0, 0.0);

    //Initialisation aléatoire de la position de la chaîne
    size_t site_current = random_site(rng);
    int mu_current = random_dir(rng);
    int epsilon_current = 2 * random_eps(rng) -1;

    //Initialisation aléatoire des theta limites pour sample et refresh
    double theta_sample{};
    double theta_refresh{};
    if (poisson) {
        theta_sample = random_theta_sample(rng);
        theta_refresh = random_theta_refresh(rng);
    }
    else {
        theta_sample = param_theta_sample;
        theta_refresh = param_theta_refresh;
    }

    //Initialisation des angles totaux parcourus à 0.0
    double theta_parcouru_sample = 0.0;
    double theta_parcouru_refresh= 0.0;

    //Angle d'update
    double theta_update = 0.0;

    //Arrays utilisés à chaque étape de la chaîne (évite de les initialiser des milliers de fois)
    array<double,6> reject_angles = {0.0, 0.0, 0.0, 0.0, 0.0};
    array<SU3,6> list_staple;

    SU3 R = random_su3(rng);
    cout << "beta = " << beta << endl;

    int samples = 0;
    array<double,2> deltas = {0.0,0.0};
    size_t event_counter = 0;
    vector<double> meas_plaquette;
    //TODO : vérifier ce qu'il se passe autour de theta_sample/theta_refresh -> effets de bord
    while (samples < N_samples) {
        compute_list_staples(links, lat, site_current, mu_current, list_staple);
        compute_reject_angles(links, site_current, mu_current, list_staple, R, epsilon_current,beta,reject_angles,rng);
        auto it = std::min_element(reject_angles.begin(), reject_angles.end());
        auto j = distance(reject_angles.begin(), it); //theta_reject = reject_angles[j]
        //cout << "Angle reject : " << reject_angles[j] << endl;
        deltas[0] = theta_sample - reject_angles[j] - theta_parcouru_sample;
        deltas[1] = theta_refresh - reject_angles[j] - theta_parcouru_refresh;

        auto it_deltas = std::min_element(deltas.begin(), deltas.end());
        auto F = distance(deltas.begin(), it_deltas);

        if ((deltas[0]<0)&&(deltas[1]<0)) {
            if (F == 0) {
                //On sample
                theta_update = theta_sample - theta_parcouru_sample;
                ecmc_update(links, site_current, mu_current, theta_update, epsilon_current, R);
                cout << "Sample " << samples << ", ";
                auto plaq = plaquette_stats(links, lat);
                cout << "<P> = " << plaq.mean << " +- " << plaq.stddev << ", " << event_counter << " events" << endl;
                event_counter = 0;
                meas_plaquette.emplace_back(plaq.mean);
                samples++;
                theta_parcouru_sample = 0;
                if (poisson) theta_sample = random_theta_sample(rng); //On retire un nouveau theta_sample
                theta_parcouru_refresh += theta_update;
                //On update jusqu'au refresh
                theta_update = theta_refresh - theta_parcouru_refresh;
                ecmc_update(links, site_current, mu_current, theta_update, epsilon_current, R);
                theta_parcouru_sample += theta_update;
                theta_parcouru_refresh = 0;
                if (poisson) theta_refresh = random_theta_refresh(rng); //On retire un nouveau theta refresh
                //On refresh
                event_counter++;
                site_current = random_site(rng);
                mu_current = random_dir(rng);
                epsilon_current = 2* random_eps(rng) -1;
                R = random_su3(rng);
            }
            if (F == 1) {
                //On update jusqu'au refresh
                theta_update = theta_refresh - theta_parcouru_refresh;
                ecmc_update(links, site_current, mu_current, theta_update, epsilon_current, R);
                theta_parcouru_sample += theta_update;
                theta_parcouru_refresh = 0;
                if (poisson) theta_refresh = random_theta_refresh(rng); //On retire un nouveau theta_refresh
                //On refresh
                event_counter++;
                site_current = random_site(rng);
                mu_current = random_dir(rng);
                epsilon_current = 2* random_eps(rng) -1;
                R = random_su3(rng);
            }
        }
        else if (deltas[F]<0) {
            if (F == 0) {
                //On update jusqu'a theta_sample
                theta_update = theta_sample - theta_parcouru_sample;
                ecmc_update(links, site_current, mu_current, theta_update, epsilon_current, R);
                //On sample
                cout << "Sample " << samples << ", ";
                auto plaq = plaquette_stats(links, lat);
                cout << "<P> = " << plaq.mean << " +- " << plaq.stddev << ", " << event_counter << " events" << endl;
                event_counter = 0;
                meas_plaquette.emplace_back(plaq.mean);
                samples++;
                theta_parcouru_sample = 0;
                if (poisson) theta_sample = random_theta_sample(rng); //On retire un nouveau theta_sample
                theta_parcouru_refresh += theta_update;
                //On finit l'update et on lift
                theta_update = -deltas[F];
                ecmc_update(links, site_current, mu_current, theta_update, epsilon_current, R);
                theta_parcouru_sample += theta_update;
                theta_parcouru_refresh += theta_update;
                //On lifte
                event_counter++;
                auto l = lift(links, lat, site_current, mu_current, j, R, lambda_3, rng);
                site_current = l.first.first;
                mu_current = l.first.second;
                epsilon_current = l.second;
            }
            if (F==1) {
                //On update jusqu'à theta_refresh
                theta_update = theta_refresh - theta_parcouru_refresh;
                ecmc_update(links, site_current, mu_current, theta_update, epsilon_current, R);
                theta_parcouru_sample += theta_update;
                theta_parcouru_refresh = 0;
                if (poisson) theta_refresh = random_theta_refresh(rng); //On retire un nouveau theta_refresh
                //On refresh
                event_counter++;
                site_current = random_site(rng);
                mu_current = random_dir(rng);
                epsilon_current = 2* random_eps(rng) -1;
                R = random_su3(rng);
            }
        }
        else {
            //On update
            theta_update = reject_angles[j];
            ecmc_update(links, site_current, mu_current, theta_update, epsilon_current, R);
            theta_parcouru_sample += theta_update;
            theta_parcouru_refresh += theta_update;
            //On lift
            event_counter++;
            auto l = lift(links, lat, site_current, mu_current, j, R, lambda_3, rng);
            site_current = l.first.first;
            mu_current = l.first.second;
            epsilon_current = l.second;
        }
    }
    return meas_plaquette;
}

int update_until_reject_d(vector<Complex> &links, size_t site, int mu, array<SU3, 6> &list_staple,
    const SU3 &R, int epsilon, const double &beta, const double &eta, mt19937_64 &rng) {
    //Propose au lien des updates avec embedding jusqu'à ce qu'exactement une plaquette refuse le move.
    //Continue les propositions tant que plusieurs ou aucune plaquette refuse.
    //Renvoie l'angle de rejet et l'indice de la plaquette dans list_staples correspondant
    double theta_update = 0.0;
    array<int, 6> accepted_angles;
    uniform_real_distribution<double> theta_dist(0.2, eta); //Pour tirer l'angle theta entre 0 et eta
    uniform_real_distribution<double> unif(0.0, 1.0); //Pour test metropolis
    SU3 proposition;
    SU3 Uold;
    SU3 Unew;
    int proposed = 0;
    while (true) {
        theta_update = epsilon * theta_dist(rng);
        for (int i = 0; i < 6; i++) {
            proposition = R * el_3(theta_update) * R.adjoint();
            Uold = view_link_const(links, site, mu);
            Unew = proposition * view_link_const(links, site, mu);
            double old_tr = (Uold * list_staple[i]).trace().real();
            double new_tr = (Unew * list_staple[i]).trace().real();
            double dS = -(beta/3.0) * (new_tr - old_tr);
            proposed++;
            bool accept = false;
            if ((dS <= 0.0)||(unif(rng) < exp(-dS))) accept = true;
            accepted_angles[i] = accept ? 1 : 0;

        }
        int sum = 0;
        for (int i = 0; i < 6; i++) {
            sum += accepted_angles[i];
        }
        if (sum == 6) {
            ecmc_update(links, site, mu, theta_update, epsilon, R);
        }
        if (sum == 5) {
            //cout << proposed << " propositions" << endl;
            auto it  = min_element(accepted_angles.begin(), accepted_angles.end());
            int j = distance(accepted_angles.begin(), it);
            return j;
        }
    }
}
