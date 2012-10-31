#ifndef CHARODE_DDS_NEW_H
#define CHARODE_DDS_NEW_H

#include <misc.h>
#include <grid.h>

class charODE_dds_new
{
    public:
        charODE_dds_new(GridClass& grid, const double mu);
        void operator() (const std::vector<double>& r,
                         std::vector<double>&       dr_ds,
                         const double          s);

        friend double gamma_ltz(double beta);

    private:
        GridClass* m_grid;
        double m_mu;
};

#endif
