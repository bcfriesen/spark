#ifndef CHARODE_DDR_H
#define CHARODE_DDR_H

#include <misc.h>
#include <grid.h>

class charODE_ddr
{
    public:
        charODE_ddr(GridClass& grid);
        void operator() (const std::vector<double>& x,
                         std::vector<double>&       dxdr,
                         const double               r);

        friend double gamma_ltz(double beta);

    private:
        GridClass* m_grid;
};

#endif
