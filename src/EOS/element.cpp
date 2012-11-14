#include <element.h>
#include <my_exceptions.h>

ElementClass::ElementClass(int Z, int I)
{
    if (I > Z + 1) throw InvalidElement(Z, I);
    /* Calculate Kurucz ion ID from Z and I. */
    const int kurucz_id = Z * 100 + I;

    // TODO: calculate partition function
}
