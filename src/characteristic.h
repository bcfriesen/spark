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
        /** Get index of tangent layer. */
        int get_tangent_layer_index();

        void push_s(double s_);
        void push_mu(double mu_);
        double get_s(int i);
        double get_mu(int i);

        friend double gamma_ltz(double beta);

    private:
        /** Impact parameter. */
        double m_p;
        /** Path length along ray (function of \f$r\f$). */
        std::vector<double> s;
        /** Direction cosine along ray (function of \f$r\f$). */
        std::vector<double> mu;
        /** Grid. */
        GridClass* m_grid;
        /** Index of tangent layer. */
        int tangent_layer_index;
};

#endif
