#ifndef CALC_RAYS
#define CALC_RAYS

#include <grid.h>
#include <characteristic.h>
#include <misc.h>

/** \brief Integrate characteristic ray ODEs.
 *
 * Calculate \f$s(r)\f$ and \f$mu(r)\f$ by integrating the \f$d/ds\f$ and
 * \f$d/dr\f$ ODEs derived in Mihalas (1980). */
void calc_rays(GridClass &grid, std::vector<Characteristic>::iterator it_char, RayIntDir direction);

#endif
