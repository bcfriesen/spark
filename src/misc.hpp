#ifndef MISC_H
#define MISC_H

#include <vector>
#include <string>
#include <cstdlib>

/** Calculates Lorentz factor \f$ \gamma \equiv (1 - \beta^2)^{-1/2} \f$. */
double gamma_ltz(double beta);

/** Linear interpolator. I stole this whole routine from Daniel Fleischman on
 * StackOverflow. */
double interpolate(std::vector< std::pair <double, double> > table, double x);

/** Retrieves environment variables. */
std::string getEnvVar(std::string const &key);

/** Generic Gaussian function of the form \f$ f(x) \equiv a
 * e^{\frac{-(x-b)^2}{2c^2}} \f$. We'll use it to calculate intrinsic line
 * profiles. */
double gaussian(double x, double a, double b, double c);

#endif
