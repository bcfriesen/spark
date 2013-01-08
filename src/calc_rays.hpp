#ifndef CALC_RAYS_H
#define CALC_RAYS_H

#include <grid.hpp>
#include <characteristic.hpp>
#include <boost/numeric/odeint.hpp>

/** \brief Integrate "front" characteristic ray ODEs.
 *
 * Calculate \f$ s(r) \f$ and \f$ mu(r) \f$ by integrating the \f$ d/dr \f$ ODE
 * derived in Hauschildt (1992).
 * */
void calc_rays(GridClass* grid, std::vector<CharNCI_F> &ray_vector);

/** \brief Integrate "back" characteristic ray ODEs.
 *
 * Calculate \f$ s(r) \f$ and \f$ mu(r) \f$ by integrating the \f$ d/dr \f$ ODE
 * derived in Hauschildt (1992).
 * */
void calc_rays(GridClass* grid, std::vector<CharNCI_B> &ray_vector);

#endif // CALC_RAYS_H
