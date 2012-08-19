#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <vector>

using namespace std;

//! Interpolator
/**
 * Right now it does only linear interpolation.
 */
double interpolate(vector<pair <double, double> > table, double x);

#endif
