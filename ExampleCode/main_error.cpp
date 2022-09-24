// Example code for the meeting of the UniPD cosmology group, Asiago 2022.
// Calculates the comoving distance to a set redshift in flat LCDM without radiation.

#include <iostream>

#include "background.h"

const double c_in_km_s = 2.99792458e5;

#include <boost/math/quadrature/gauss_kronrod.hpp>
using boost::math::quadrature::gauss_kronrod;



int main()
{
// Set parameters
    double z = 1000.;
    double z0 = 0.;
    
    ACDM::Background background;
    background.Om0 = 0.3;
    background.H0_in_km_over_s_Mpc = 70.;
    
// Define the integrand
    auto comoving_distance_integrand = [&] (const double z)
    {
        return c_in_km_s / sqrt(background.E2(z)) / background.H0_in_km_over_s_Mpc;
    };
    
// Integrate comoving distance
    double error;
    const double epsrel = 1.e-6;
    const int max_depth = 1;
    double chi = gauss_kronrod<double, 5>::integrate(comoving_distance_integrand,
                                                      z0, z,
                                                      max_depth, epsrel, &error);                     
// Remove error management because why not    
//     if(error/chi > epsrel)
//     {
//         std::cerr << "Integration failed with I = " << chi << ", error = " << error << std::endl;
//         return 1;
//     }
    
    
// Output
    std::cout << "The comoving distance at z=" << z << " is chi=" << chi << " Mpc (correct answer is ~14 Gpc)\n";
    
    return 0;
}
