#ifndef CHARODE_DDS_H
#define CHARODE_DDS_H

#include <misc.h>
#include <grid.h>

/** \brief Characteristic ray \f$d/ds\f$ ODEs from Mihalas (1980).
 *
 * ODE class supplied to integrator to solve \f$r(s)\f$ and \f$\mu(s)\f$.
 */
class charODE_dds
{
    public:
        /** Use grid data to set up parameters. */
        charODE_dds(GridClass& grid);
        /** This is the function called by odeint. */
        void operator() (const std::vector<double>& x,
                         std::vector<double>&       dxds,
                         const double               s);

        friend double gamma_ltz(double beta);

    private:
        GridClass* m_grid;
};

#endif
