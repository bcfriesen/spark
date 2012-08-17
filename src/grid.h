#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <vector>

using namespace std;

class GridClass
{
    public:
        GridClass();
        ~GridClass();
        void initialize(string filename);
        double rad(int layer) const; // read from tabulated data in layer file
        double vel(int layer) const; // read from tabulated data in layer file
        double vel(double rad) const; // interpolate
        double beta(int layer) const; // read from tabulated data in layer file
        double beta(double rad) const; // interpolate
        double dbeta_dr(int layer) const; // read from tabulated data in layer file
        double dbeta_dr(double rad) const; // interpolate
        friend double interpolate(vector< pair<double, double> > table, double x);
        friend double gamma_ltz(double beta);
    private:
        vector< pair<double, double> > rad_vel; // radius and velocity of each layer
};

#endif
