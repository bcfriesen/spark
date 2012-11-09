#ifndef CALC_RAYS
#define CALC_RAYS

#include <grid.h>
#include <characteristic.h>
#include <misc.h>

void calc_rays(GridClass &grid, std::vector<Characteristic>::iterator it_char, RayIntDir direction);

#endif
