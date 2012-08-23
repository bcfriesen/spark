#ifndef MISC_H
#define MISC_H

//! Calculates Lorentz factor \f$ \gamma \equiv (1 - \beta^2)^{-1/2} \f$
double gamma_ltz(double beta)
{
    return 1.0 / sqrt(1.0 - pow(beta, 2));
}

#endif
