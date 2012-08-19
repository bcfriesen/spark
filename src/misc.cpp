#ifndef MISC_H
#define MISC_H

#include <math.h>
#include <misc.h>

//! Calculates \f$ \beta \equiv v/c \f$
double gamma_ltz(double beta)
{
    return 1.0 / sqrt(1.0 - pow(beta, 2));
}

#endif
