#ifndef COMBINE_TANGENT_RAYS_HPP
#define COMBINE_TANGENT_RAYS_HPP

#include <characteristic.hpp>

/** \brief Combine "front" and "back" tangent rays to form a single ray.
 *
 * The ODEs for calculating co-moving quantities along the front and back rays
 * look slightly different, thanks to a few minus signs, which is why we
 * calculate the front and back half separately.  But ultimately we want to
 * treat each tangent ray as a single entity all the way through the
 * atmosphere, so after we calculate the front and back half of each ray, we
 * will send them here to be stitched together at \f$ s = 0 \f$, the center of
 * the atmosphere.
 * */

CharNCI combine_tangent_rays(const CharNCI_B back_ray, const CharNCI_F front_ray, GridClass& grid);

#endif // COMBINE_TANGENT_RAYS_HPP
