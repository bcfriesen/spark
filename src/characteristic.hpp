#ifndef CHARACTERISTIC_H
#define CHARACTERISTIC_H

#include <vector>
#include <grid.hpp>

/** \brief Characteristic ray class.
 *
 * Holds values of \f$ s(r) \f$ and \f$ \mu(r) \f$ for each ray calculated by
 * using Bin Chen's RT formalism in Chen et al. (2007).
 * */

class Characteristic
{
    public:
        /** Constructor requires a grid, otherwise our characteristics have no structure. */
        Characteristic(GridClass& grid);
        /** Save resulting values of \f$ s \f$ and \f$ \mu \f$ as we integrate
         * from one radial point to the next. */
        void push_s_mu(double s, double mu);
        /** Get the geometric coordinate \f$ s \f$ at point \f$ k \f$ along ray
         * \f$ i \f$. */
        double get_s(int k) const;
        /** Get the direction cosine \f$ \mu \f$ at point \f$ k \f$ along ray
         * \f$ i \f$. */
        double get_mu(int k) const;
        /** Return total number of points along ray. */
        unsigned int get_num_ray_pts() const;
        /** Iterator lower bound for the vector of \f$ s \f$ points. */
        std::vector<double>::const_iterator s_vec_begin() const;
        /** Iterator upper bound for the vector of \f$ s \f$ points. */
        std::vector<double>::const_iterator s_vec_end() const;
        /** Iterator lower bound for the vector of \f$ mu \f$ points. */
        std::vector<double>::const_iterator mu_vec_begin() const;
        /** Iterator upper bound for the vector of \f$ mu \f$ points. */
        std::vector<double>::const_iterator mu_vec_end() const;

    protected:
        /** Path length along ray. */
        std::vector<double> m_s;
        /** Direction cosine along ray. */
        std::vector<double> m_mu;
        /** Pointer to grid variables. */
        GridClass* m_grid;

    friend double gamma_ltz(double beta);
};

/** Non-core-intersecting characteristic ray. */
class TangentRay: public Characteristic
{
    public:
        /** Initialize a ray tangent to layer \f$ j \f$. For simplicity we
         * define the geometric coordinate at this point as \f$ s=0 \f$. */
        TangentRay(GridClass& grid, int i);
        /** Return the impact parameter for this characteristic ray at \f$ s=0
         * \f$. This is the same as just the radius of the tangent layer. */
        double get_p() const;
        /** Return index of tangent layer. */
        unsigned int tangent_layer_index() const;

    protected:
        /** Impact parameter at \f$ s=0 \f$. This is the same as just the
         * radius of the tangent layer. */
        double m_p;
        /** Index of tangent layer. */
        unsigned int m_tangent_layer_index;
};

// TODO: What else do core-intersecting rays need?
/** Core-intersecting characteristic ray. */
class CoreIntersectingRay: public Characteristic
{
    public:
        /** Initialize a ray emerging from the core at impact parameter \f$ p \f$. Can't user
         * layer indices anymore because we're below the first layer by definition. */
        CoreIntersectingRay(GridClass& grid, double p);
};

#endif
