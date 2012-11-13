#ifndef CHARACTERISTIC_H
#define CHARACTERISTIC_H

#include <vector>
#include <grid.h>

/** \brief Characteristic ray class.
 *
 * Holds values of \f$s(r)\f$ and \f$\mu(r)\f$ for each ray calculated by
 * integrating the characteristic ray ODEs in Mihalas (1980). */
class Characteristic
{
    public:
        /** Initialize a ray tangent to layer \f$i\f$ at \f$s=0\f$. */
        Characteristic(GridClass& grid, int i);
        /** Return the impact parameter for this characteristic ray at \f$s=0\f$. */
        double get_p();
        /** Return index of tangent layer. */
        int get_tangent_layer_index();
        /** Save the resulting value of \f$s\f$ as we integrate from one radial
         * point to the next. */
        void push_s(double s_);
        /** Save the resulting value of \f$\mu\f$ as we integrate from one radial
         * point to the next. */
        void push_mu(double mu_);
        /** Return \f$s(r_i)\f$. */
        double get_s(int i);
        /** Return \f$\mu(r_i)\f$. */
        double get_mu(int i);

        friend double gamma_ltz(double beta);

    private:
        /** Impact parameter at \f$s=0\f$. */
        double m_p;
        /** Path length along ray (function of \f$r\f$). */
        std::vector<double> s;
        /** Direction cosine along ray (function of \f$r\f$). */
        std::vector<double> mu;
        /** Pointer to grid variables. */
        GridClass* m_grid;
        /** Index of tangent layer. */
        int tangent_layer_index;
};

#endif
