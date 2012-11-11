#ifndef GRID_H
#define GRID_H

#include <vector>

/** Contains most grid-based variables and functions. */
class GridClass
{
    public:
        /** Reads the radius/velocity grid data from the layer file and stores
         *  it in member variables. */
        GridClass(char* yaml_file);
        /** Get number of layers in layer file. */
        int get_num_layers() const;
        /** Get number of core-intersecting rays (this is a knob). */
        int get_num_core_intersect_rays() const;
        /** Get radius from tabulated data in layer file. */
        double rad(int layer) const;
        /** Get velocity from tabulated data in layer file. */
        double vel(int layer) const;
        /** Interpolate velocity between tabulated data points. */
        double vel(double rad) const;
        /**  Get \f$ v/c \f$ from tabulated data in layer file. */
        double beta(int layer) const;
        /** Interpolate \f$ v/c \f$ from tabulated data in layer file. */
        double beta(double rad) const;
        /** Calculate \f$ d\beta/dr \f$ at tabulated data points. */
        double dbeta_dr(int layer) const;
        /** Calculate \f$ d\beta/dr \f$ by interpolating between data points. */
        // FIXME: this whole thing is probably completely wrong
        double dbeta_dr(double rad) const;

        /* Sometimes it's useful to iterate over the grid variables from
         * outside this class. */
        std::vector< std::pair<double, double> >::iterator begin();
        std::vector< std::pair<double, double> >::iterator end();

        friend double interpolate(std::vector< std::pair<double, double> > table, double x);
        friend double gamma_ltz(double beta);

    private:
        /** Radius and velocity of each layer. */
        std::vector< std::pair<double, double> > rad_vel;
        /** Number of core-intersecting rays. */
        int num_core_intersect_rays;
};

#endif
