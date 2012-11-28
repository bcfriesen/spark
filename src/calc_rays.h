#ifndef CALC_RAYS_H
#define CALC_RAYS_H

#include <grid.h>
#include <characteristic.h>
#include <boost/numeric/odeint.hpp>

// Functions with template parameters must be defined in the header file - you
// can't split them up.

template <typename T>
/** \brief Integrate characteristic ray ODEs.
 *
 * Calculate \f$ s(r) \f$ and \f$ mu(r) \f$ by integrating the \f$ d/dr \f$ ODE
 * derived in Hauschildt (1992).
 * */
void calc_rays(const GridClass &grid, std::vector<T> &ray_vector)
{
    // Iterates over vector of characteristic rays.
    typename std::vector<T>::iterator it_char;

    /* s coordinate for ds/dr equation. Even for a single ODE, we must supply
     * the variable as a vector because odeint uses iterators to keep track of
     * variables. */
    std::vector<double> s_of_r(1);

    // I read somewhere on StackOverflow that it's faster to use ++blah instead
    // of blah++, especially for iterators which track vectors of complicated
    // classes, because the latter requires a call to the copy constructor
    // whereas the former does not.
    for (it_char = ray_vector.begin(); it_char != ray_vector.end(); ++it_char)
    {
        /* Integrate ds/dr a la Hauschildt (1992). */
        /* By definition s = 0 at the tangent point. */
        s_of_r.at(0) = 1.0e-30*grid.rad(0);
        for (int i = it_char->get_tangent_layer_index(); i < grid.get_num_layers()-1; i++)
        {
            if (i == it_char->get_tangent_layer_index())
            {
                boost::numeric::odeint::integrate(*it_char, s_of_r, (1.0 + 1.0e-10)*grid.rad(i), grid.rad(i+1), 1.0e-5*grid.rad(i));
            }
            else
            {
                boost::numeric::odeint::integrate(*it_char, s_of_r, grid.rad(i), grid.rad(i+1), 1.0e-5*grid.rad(i));
            }
            // TODO: fix the indices on s(r_i) in the characteristic class so
            // they match the indices on the vector of radial points {r_i} in
            // the grid class.
            // TODO: figure out how to export the analytic values of mu(r_i)
            // calculated in the overloaded () operator to the member vector
            // of the characteristic class. The way it is now, it's calculated
            // twice: once in the () operator, and again here.
            const double mu_E = it_char->sign_of_mu() * sqrt(1.0 - (pow(it_char->get_p(), 2) / pow(grid.rad(i+1), 2)));
            const double mu   = (mu_E - grid.beta(i+1)) / (1.0 - grid.beta(i+1) * mu_E);
            it_char->push_s_mu(s_of_r.at(0), mu);
        }
    }
}

#endif
