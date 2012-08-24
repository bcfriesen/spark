#ifndef CHARACTERISTIC_H
#define CHARACTERISTIC_H

#include <grid.h>

class CharacteristicODEClass
{
    public:
        CharacteristicODEClass(GridClass &grid, int i);
        void operator() (const std::vector<double> &x,
                         std::vector<double>       &dxds,
                         const double         s);
        //! get impact parameter for this ray
        double get_p();

        friend double gamma_ltz(double beta);

    private:
        double m_p;/*!< impact parameter */
        GridClass* m_grid; /*!< grid */
};

#endif
