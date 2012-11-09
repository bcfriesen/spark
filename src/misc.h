#ifndef MISC_H
#define MISC_H

#include <vector>

//! Calculates Lorentz factor \f$ \gamma \equiv (1 - \beta^2)^{-1/2} \f$
double gamma_ltz(double beta);
/**
 * Linear interpolator.
 * I stole this whole routine from Daniel Fleischman on StackOverflow.
 */
double interpolate(std::vector<std::pair <double, double> > table, double x);

//! Which way do we integrate the characteristic rays?
enum RayIntDir {FORWARD = 0, BACKWARD = 1};

#endif
