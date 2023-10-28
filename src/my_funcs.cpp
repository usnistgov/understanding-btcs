#include <math.h>  // todo: do I need this?
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <stdio.h>
#include <iostream>
//#include <boost/math/quadrature>

using namespace std;

namespace bmath = boost::math;
//using namespace boost::math::quadrature;
//using boost::multiprecision::cpp_bin_float_quad;

extern "C" double exp_I0(double T, double z){
    return exp(-z - T) * bmath::cyl_bessel_i(0, 2 * sqrt(z*T));
}

extern "C" double integrate(double T_jm1, double T_j, double z){
    size_t max_refinements = 15;

    double I = bmath::quadrature::trapezoidal(
        std::bind(exp_I0, std::placeholders::_1, z), T_jm1, T_j, 1.e-9, max_refinements
        );
    return I;
}

extern "C" void update_solution(
    double *c_j, double *c_jm1, double T_jm1, double T_j, double *X, int n
){
    int i = 0;
    c_j[0] = 1;
    for (i = 1; i < n; i++)
    {
        c_j[i] = c_jm1[i] \
                    + exp_I0(T_j, X[i]) \
                    - exp_I0(T_jm1, X[i]) \
                    + integrate(T_jm1, T_j, X[i]);
    }
}


extern "C" void analytical_btc(
    double *c, double *T, double X, int m
){
    int j = 0;
    while (T[j] <= 0.)
    {
        c[j++] = 0;
    }

    c[j++] = exp_I0(T[j], X) + integrate(0, T[j], X);
    while (j < m)
    {
        c[j++] = c[j-1] \
                    + exp_I0(T[j], X) \
                    - exp_I0(T[j-1], X) \
                    + integrate(T[j-1], T[j], X);
    }
}