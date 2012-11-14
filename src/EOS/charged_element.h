#ifndef CHARGED_ELEMENT_H
#define CHARGED_ELEMENT_H

#include <element.h>

class ChargedElementClass: public ElementClass
{
    public:
        /** Chargd elements can have any ionization state \f$ I = 0 \f$. */
        ChargedElementClass(int Z, int I);

    private:
};

#endif
