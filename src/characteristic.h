#ifndef CHARACTERISTIC_H
#define CHARACTERISTIC_H

#include <vector>
#include <grid.h>

/** Characteristic ray. */
class Characteristic
{
    public:
        Characteristic(GridClass& grid, int i);
        /** Return the impact parameter for this characteristic ray. */
        double get_p();
        void push_s(double s_);
        void push_mu(double mu_);
        friend double gamma_ltz(double beta);
        double get_s(int i);
        double get_mu(int i);

    private:
        double m_p; /** Impact parameter. */
        std::vector<double> s; /** Path length along ray (function of \f$r\f$). */
        std::vector<double> mu; /** Direction cosine along ray (function of \f$r\f$). */
        GridClass* m_grid; /** Grid. */
};

#endif
