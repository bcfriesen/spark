#ifndef NEUTRAL_ELEMENT_H
#define NEUTRAL_ELEMENT_H

class NeutralElementClass: public ElementClass
{
    public:
        /** Neutral elements have ionization state \f$ I = 0 \f$. */
        NeutralElementClass(int Z);

    private:
};

#endif
