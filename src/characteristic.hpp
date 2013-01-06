#ifndef CHARACTERISTIC_H
#define CHARACTERISTIC_H

#include <vector>
#include <grid.hpp>

/** \brief Characteristic ray class.
 *
 * Holds values of \f$ s(r) \f$ and \f$ \mu(r) \f$ for each ray calculated by
 * integrating the characteristic ray ODEs in Mihalas (1980).
 * */

class Characteristic
{
    public:
        /** Constructor requires a grid, otherwise our characteristics have no structure. */
        Characteristic(GridClass& grid);
        /** Save resulting values of \f$ s \f$ and \f$ \mu \f$ as we integrate
         * from one radial point to the next. */
        void push_s_mu(double s, double mu);
        /** Return \f$ s(r_i) \f$. */
        double get_s(int i);
        /** Return \f$ \mu(r_i) \f$. */
        double get_mu(int i);
        /** Forward and backward rays are distinguished by the sign of \f$ mu(r)
         * \f$ (see Eqs. 14-15 of Hauschildt (1992)). */
        virtual double sign_of_mu() = 0;
        /** Iterator lower bound for the vector of \f$ s(r_i) \f$ and \f$
         * mu(r_i) \f$ pairs. */
        std::vector< std::pair<double, double> >::const_iterator s_mu_vec_begin();
        /** Iterator upper bound for the \f$ s(r_i) \f$ vector. */
        std::vector< std::pair<double, double> >::const_iterator s_mu_vec_end();

    protected:
        /** Path length and direction cosine along ray (both functions of \f$ r
         * \f$). */
        std::vector< std::pair<double, double> > m_s_mu;
        /** Pointer to grid variables. */
        GridClass* m_grid;

    friend double gamma_ltz(double beta);
};

/** Non-core-intersecting characteristic ray. */
class CharNCI: public Characteristic
{
    public:
        /** Initialize a ray tangent to layer \f$ i \f$ at \f$ s=0 \f$. */
        CharNCI(GridClass& grid, int i);
        /** Return the impact parameter for this characteristic ray at \f$ s=0 \f$. */
        double get_p();
        /** Return index of tangent layer. */
        int get_tangent_layer_index();

    protected:
        /** Impact parameter at \f$ s=0 \f$. */
        double m_p;
        /** Index of tangent layer. */
        int m_tangent_layer_index;
};

/** Non-core-intersecting rays on the "front" side of the grid, closer to the
 * observer. */
class CharNCI_F: public CharNCI
{
    public:
        CharNCI_F(GridClass& grid, int i);
        double sign_of_mu();
        /** This is the function called by odeint. */
        void operator() (const std::vector<double>& x,
                         std::vector<double>&       dsdr,
                         const double               r);
};

/** Non-core-intersecting rays on the "back" side of the grid, farther from the
 * observer. */
class CharNCI_B: public CharNCI
{
    public:
        CharNCI_B(GridClass& grid, int i);
        double sign_of_mu();
        /** This is the function called by odeint. */
        void operator() (const std::vector<double>& x,
                         std::vector<double>&       dsdr,
                         const double               r);
};

// TODO: What else do core-intersecting rays need?
/** Core-intersecting characteristic ray. */
class CharCI: public Characteristic
{
    public:
        /** Initialize a ray emerging from the core at impact parameter \f$ p \f$. Can't user
         * layer indices anymore because we're below the first layer by definition. */
        CharCI(GridClass& grid, double p);
};

/** Core-intersecting rays on the "front" side of the grid, closer to the
 * observer. */
class CharCI_F: public CharCI
{
    public:
        CharCI_F(GridClass& grid, double p);
        double sign_of_mu();
        /** This is the function called by odeint. */
        void operator() (const std::vector<double>& x,
                         std::vector<double>&       dsdr,
                         const double               r);
};

/** Core-intersecting rays on the "back" side of the grid, farther from the
 * observer. */
class CharCI_B: public CharCI
{
    public:
        CharCI_B(GridClass& grid, double p);
        double sign_of_mu();
        /** This is the function called by odeint. */
        void operator() (const std::vector<double>& x,
                         std::vector<double>&       dsdr,
                         const double               r);
};

#endif
