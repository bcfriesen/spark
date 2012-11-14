#ifndef ELEMENT_H
#define ELEMENT_H

#include <species.h>

/** \brief A non-electron species. */
class ElementClass: public SpeciesClass
{
    public:
        /** All elements (besides electrons) have an atomic number and an
         * ionization stage. The nomenclature here follows Hummer & Mihalas
         * (1988): \f$ I = 1 \f$ is neutral. */
        ElementClass(int Z, int I);

    protected:
        /** Atomic number. */
        int m_Z;
        /** Ionization state (\f$ I = 1 \f$ is neutral). */
        int m_I;
};

#endif
