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
    public:
        /** Set number abundance. */
        void set_abund(double N);
        /** Return number abundance. */
        double get_abund() const;

    protected:
        /** Particle mass. */
        /* TODO: Implement function to set particle mass. This will come from
         * reading model atom data. But first I need to learn HDF5. */
        double m_M_s;
        /** Number abundance. */
        double m_N;
        /** Partition function. */
        double m_part_Z;
};

#endif
