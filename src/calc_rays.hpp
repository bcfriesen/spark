#ifndef CALC_RAYS_H
#define CALC_RAYS_H

#include <grid.hpp>
#include <characteristic.hpp>

/** \brief Calculate geometric properties of characteristic rays.
 *
 * The beauty of Bin's RT formalism is that the geometric properties of the
 * characteristic rays are calculated in the *observer's* frame, which means
 * they're straight lines in flat spacetime. So instead of calculating \f$ s(r)
 * \f$ and \f$ mu(r) \f$ by integrating the \f$ d/dr \f$ ODE derived in
 * Hauschildt (1992) and/or Mihalas (1980), both of which use coordinates in
 * mixed frames to track the characteristics, we just draw some triangles, use
 * high some high school geometry, and we're done! Awesome!
 * */
void calc_rays(GridClass* grid, std::vector<TangentRay> &ray_vector);

#endif // CALC_RAYS_H
