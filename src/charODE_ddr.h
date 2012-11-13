#ifndef CHARODE_DDR_H
#define CHARODE_DDR_H

#include <misc.h>
#include <grid.h>

/** \brief Characteristic ray \f$d/dr\f$ ODEs from Mihalas (1980).
 *
 * ODE class supplied to integrator to solve \f$s(r)\f$ and \f$\mu(r)\f$.
 */
class charODE_ddr
{
    public:
        /** Use grid data to set up parameters. */
        charODE_ddr(GridClass& grid);
        /** This is the function called by odeint. */
        void operator() (const std::vector<double>& x,
                         std::vector<double>&       dxdr,
                         const double               r);

        friend double gamma_ltz(double beta);

    private:
        /** Pointer to grid data. */
        GridClass* m_grid;
};

#endif
