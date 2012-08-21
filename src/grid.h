#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <vector>

using namespace std;

/**
 * Contains most grid-based variables and functions. Which is pretty much
 * everything.
 */
class GridClass
{
    public:
        //! Does most initialization.
        /**
         * Reads the radius/velocity grid data from the layer file and stores it
         * in member variables.
         */
        GridClass(char* layer_file);
        ~GridClass();
        double rad(int layer) const; //!< Get radius from tabulated data in layer file.
        double vel(int layer) const; //!< Get velocity from tabulated data in layer file.
        double vel(double rad) const; //!< Interpolate velocity between tabulated data points.
        double beta(int layer) const; //!<  Get \f$ v/c \f$ from tabulated data in layer file.
        double beta(double rad) const; //!< Interpolate \f$ v/c \f$ from tabulated data in layer file.
        double dbeta_dr(int layer) const; //!< Calculate \f$ d\beta/dr \f$ at tabulated data points.
        double dbeta_dr(double rad) const; //!< Calculate \f$ d\beta/dr \f$ by interpolating between data points.

        friend double interpolate(vector< pair<double, double> > table, double x);
        friend double gamma_ltz(double beta);
    private:
        vector< pair<double, double> > rad_vel; /*!< radius and velocity of each layer */
};

#endif
