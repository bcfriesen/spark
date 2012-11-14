#ifndef SPECIES_H
#define SPECIES_H

#include <misc.h>

/** \brief Generic particle class.
 *
 * All particles involved in calculating the EOS are some kind of species. The
 * only thing they all have in common is the partition function.
 * */
class SpeciesClass
{
    protected:
        /** Partition function. Every species gets one of these. */
        double m_part_Z;
};

#endif
