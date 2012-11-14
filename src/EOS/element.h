#ifndef ELEMENT_H
#define ELEMENT_H

#include <species.h>

class ElementClass: public SpeciesClass
{
    public:
        /* All elements (besides electrons) have an atomic number and an
         * ionization stage. The nomenclature here follows Hummer & Mihalas
         * (1988): \f$ I = 1 \f$ is neutral. */
        ElementClass(int Z, int I);
    private:
        /** Atomic number. */
        int m_Z;
        /** Ionization state (\f$ I = 1 \f$ is neutral). */
        int m_I;
};

#endif
