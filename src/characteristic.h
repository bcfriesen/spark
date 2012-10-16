#ifndef CHARACTERISTIC_H
#define CHARACTERISTIC_H

#include <vector>
#include <grid.h>

class Characteristic
{
    public:
        Characteristic(GridClass& grid, int i);
        //! Return the impact parameter for this characteristic ray.
        double get_p();
        void push_s(double s_);
        void push_mu(double mu_);
        friend double gamma_ltz(double beta);

    private:
        double m_p; /*!< impact parameter */
        std::vector<double> s; /*!< length along ray (function of \f$r\f$)*/
        std::vector<double> mu; /*!< direction-cosine along ray (function of \f$r\f$) */
        GridClass* m_grid; /*!< grid */
};

#endif
