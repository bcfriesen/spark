#ifndef CHARODE_DDR_H
#define CHARODE_DDR_H

class charODE_ddr
{
    public:
        void operator() (const std::vector<double>& x,
                         std::vector<double>&       dxds,
                         const double               s);

        friend double gamma_ltz(double beta);
};

#endif
