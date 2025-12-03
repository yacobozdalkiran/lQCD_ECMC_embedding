//
// Created by ozdalkiran-l on 11/13/25.
//
#include "su3utils.h"
#include "../lattice/lattice.h"

SU3 random_su3(std::mt19937_64 &rng) {
    //Génère une matrice de SU3 aléatoire uniformément selon la mesure de Haar en utilisant la décomposition QR
    std::normal_distribution<double> gauss(0.0, 1.0);
    Eigen::Matrix3cd z;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            z(i,j) = Complex(gauss(rng), gauss(rng));

    // Décomposition QR
    Eigen::HouseholderQR<SU3> qr(z);
    SU3 Q = qr.householderQ();

    // Corrige la phase globale pour que det(Q) = 1
    Complex detQ = Q.determinant();
    Q /= std::pow(detQ, 1.0/3.0);

    return Q;
}

SU3 su2_quaternion_to_su3(const array<double,4> &su2, int i, int j){
    //Permet d'embedder une matrice SU2 représentation quaternionique en matrice SU3 (sous groupe d'indice i,j)
    if (i==j) cerr<<"i = j wrong embedding\n";
    SU3 X;
    int k = 3-i-j;
    X.setZero();
    X(k,k) = Complex(1.0,0.0);
    X(i,i) = Complex(su2[0],su2[3]);
    X(j,j) = Complex(su2[0],-su2[3]);
    X(i,j) = Complex(su2[2],su2[1]);
    X(j,i) = Complex(-su2[2],su2[1]);
    return X;
}

SU3 random_SU3_epsilon(double epsilon, mt19937_64 &rng) {
    //Pour générer des matrices de SU3 epsilon proches de l'identité (cf Gattringer)
    uniform_real_distribution<double> unif(-0.5,0.5);
    array<double,4> x = {0.0, 0.0, 0.0, 0.0};
    SU3 M = SU3::Identity();

    //double r0 = unif(rng);
    double r1 = unif(rng);
    double r2 = unif(rng);
    double r3 = unif(rng);
    double norm = sqrt(r1*r1 + r2*r2 + r3*r3);
    x[0] = sqrt(1-epsilon*epsilon);
    x[1] = epsilon * r1 / norm;
    x[2] = epsilon * r2 / norm;
    x[3] = epsilon * r3 / norm;
    M *= su2_quaternion_to_su3(x, 0,1);

    //r0 = unif(rng);
    r1 = unif(rng);
    r2 = unif(rng);
    r3 = unif(rng);
    norm = sqrt(r1*r1 + r2*r2 + r3*r3);
    x[0] = sqrt(1-epsilon*epsilon);
    x[1] = epsilon * r1 / norm;
    x[2] = epsilon * r2 / norm;
    x[3] = epsilon * r3 / norm;
    M *= su2_quaternion_to_su3(x, 0,2);

    //r0 = unif(rng);
    r1 = unif(rng);
    r2 = unif(rng);
    r3 = unif(rng);
    norm = sqrt(r1*r1 + r2*r2 + r3*r3);
    x[0] = sqrt(1-epsilon*epsilon);
    x[1] = epsilon * r1 / norm;
    x[2] = epsilon * r2 / norm;
    x[3] = epsilon * r3 / norm;
    M *= su2_quaternion_to_su3(x, 1,2);

    return M;
}

vector<SU3> metropolis_set(double epsilon, int size, mt19937_64 &rng) {
    //Crée un set de matrices SU(3) epsilon-proches de l'identité de taille size avec leurs adjoints
    vector<SU3> set(size+1);
    set[0] = SU3::Identity();
    for (int i = 1; i < size+1; i+=2) {
        set[i] = random_SU3_epsilon(epsilon, rng);
        set[i+1] = set[i].adjoint();
    }
    return set;
}



vector<SU3> ecmc_set(double epsilon, vector<SU3> &set, mt19937_64 &rng) {
    //Crée un set de matrices SU(3) epsilon-proches de l'identité de taille size avec leurs adjoints
    size_t size = set.size()-1;
    set[0] = SU3::Identity();
    for (int i = 1; i < size+1; i+=2) {
        set[i] = random_SU3_epsilon(epsilon, rng);
        set[i+1] = set[i].adjoint();
    }
    return set;
}

void projection_su3(vector<Complex> &links, size_t site, int mu){
    //Reprojette sur SU3 en utilisant Gram-Schmidt
    SU3 U = view_link(links, site, mu);
    Eigen::Vector3cd c0 = U.col(0);
    Eigen::Vector3cd c1 = U.col(1);
    Eigen::Vector3cd c2 = U.col(2);

    c0.normalize();
    c1 -= c0 * (c0.adjoint() * c1)(0,0);
    c1.normalize();
    c2 -= c0 * (c0.adjoint() * c2)(0,0);
    c2 -= c1 * (c1.adjoint() * c2)(0,0);
    c2.normalize();

    U.col(0) = c0;
    U.col(1) = c1;
    U.col(2) = c2;

    // Fix determinant to 1
    Complex det = U.determinant();
    U /= exp(Complex(0.0, arg(det) / 3.0));
}