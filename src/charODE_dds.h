#ifndef CHARODE_DDS_H
#define CHARODE_DDS_H

#include <misc.h>
#include <grid.h>

/** Characteristic ray d/ds ODEs. */
class charODE_dds
{
    public:
        charODE_dds(GridClass& grid);
        void operator() (const std::vector<double>& x,
                         std::vector<double>&       dxds,
                         const double               s);

        friend double gamma_ltz(double beta);

    private:
        GridClass* m_grid;
};

#endif
