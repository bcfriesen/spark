#ifndef CONST_H
#define CONST_H

#include<gsl/gsl_const_cgsm.h>

/** Namespace with physical constants. Units are CGS unless explicitly stated
 *  otherwise. */
namespace physconst
{
    const double c_light = GSL_CONST_CGSM_SPEED_OF_LIGHT; /** Speed of light. */
    const double cm2km = 1.0e-5; /** Convert cm to km. */
}

#endif
