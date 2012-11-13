#include <element.h>

ElementClass::ElementClass(int Z, int I)
{
    /* Calculate Kurucz ion ID from Z and I. */
    const int kurucz_id = Z * 100 + I;
    // TODO: calculate partition function
}
