#ifndef GRID_H
#define GRID_H

#include <vector>
#include <string>

/** \brief Contains most grid-based variables and functions. */
class GridClass
{
    public:
        /** Reads the radius/velocity grid data from the layer file and stores
         *  it in member variables. */
        GridClass(std::string layer_file);
        /** Get number of layers in layer file. */
        unsigned int get_num_layers() const;
        /** Get number of core-intersecting rays (this is a knob). */
        int get_num_core_intersect_rays() const;
        /** Get radius from tabulated data in layer file. */
        double rad(unsigned int layer) const;
        /** Get velocity from tabulated data in layer file. */
        double vel(unsigned int layer) const;
        /** Interpolate velocity between tabulated data points. */
        double vel(double rad) const;
        /**  Get \f$ v/c \f$ from tabulated data in layer file. */
        double beta(unsigned int layer) const;
        /** Interpolate \f$ v/c \f$ from tabulated data in layer file. */
        double beta(double rad) const;
        /** Calculate \f$ d\beta/dr \f$ at tabulated data points. */
        double dbeta_dr(unsigned int layer) const;
        /** Calculate \f$ d\beta/dr \f$ by interpolating between data points. */
        // FIXME: this is probably completely wrong
        double dbeta_dr(double rad) const;
        /** Begin iterator for radius/velocity pair. */
        std::vector< std::pair<double, double> >::iterator begin();
        /** End iterator for radius/velocity pair. */
        std::vector< std::pair<double, double> >::iterator end();

        friend double interpolate(std::vector< std::pair<double, double> > table, double x);
        friend double gamma_ltz(double beta);

    private:
        /** Radius and velocity of each layer. Since we're assuming homologous
         * expansion, these are always paired together. */
        std::vector< std::pair<double, double> > rad_vel;

        /** Reversed radius/velocity pairs. The interpolator needs the ordinate
         * value to be monotonically increasing, but it's way more convenient
         * to track characteristic rays by starting at the outside of the
         * ejecta and moving inward, i.e., as the "index" of the rays
         * increases, the radius monotonically decreases. Instead of reversing
         * the vector of radius/velocity pairs every time we need to
         * interpolate, I'll just store a reversed copy here. This is a
         * TERRIBLE fix, but I haven't figured out a better way to do it. */
         std::vector< std::pair<double, double> > rad_vel_reversed;

        /* TODO: figure out a better way to deal with the
         * monotonically-increasing requirement of the interpolator. */

        /** Number of core-intersecting rays. */
        int num_core_intersect_rays; };

#endif
