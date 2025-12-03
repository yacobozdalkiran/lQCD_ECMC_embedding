//
// Created by ozdalkiran-l on 11/13/25.
//

#ifndef LQCD_ECMC_EMBEDDING_SU3UTILS_H
#define LQCD_ECMC_EMBEDDING_SU3UTILS_H

#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <random>

using Complex = std::complex<double>;
using SU3 = Eigen::Matrix3cd;

using namespace std;

inline int dsign(double x) {
    //fonction signe pour double
    if (x>0) return 1;
    if (x<0) return -1;
    return 0;
}

inline SU3 el_3(double xi) {
    /*Generates a exp(i xi lambda_3)*/
    SU3 result;
    double cxi = cos(xi);
    double sxi = sin(xi);
    result << Complex(cxi, sxi), Complex(0.0,0.0), Complex(0.0,0.0),
              Complex(0.0,0.0), Complex(cxi, -sxi), Complex(0.0,0.0),
              Complex(0.0,0.0), Complex(0.0,0.0), Complex(1.0,0.0);
    return result;
}

SU3 random_su3(std::mt19937_64 &rng);

SU3 su2_quaternion_to_su3(const array<double,4> &su2, int i, int j);

SU3 random_SU3_epsilon(double epsilon, mt19937_64 &rng);

vector<SU3> metropolis_set(double epsilon, int size, mt19937_64 &rng);

vector<SU3> ecmc_set(double epsilon, vector<SU3> &set, mt19937_64 &rng);

void projection_su3(vector<Complex> &links, size_t site, int mu);

#endif //LQCD_ECMC_EMBEDDING_SU3UTILS_H