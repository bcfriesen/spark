#ifndef CHARACTERISTIC_H
#define CHARACTERISTIC_H

#include <vector>

#include <grid.hpp>

/** \brief Characteristic ray data structure.
 *
 * Holds data calculated along rays (e.g., \f$ s, mu, \tau \f$).
 *
 * */
struct CharacteristicData {
    public:
        /** geometric length - defined as \f$ s=0 \f$ at the back of the ray */
        double s;
        /** direction cosine (measured in observer's frame) */
        double mu;
        /** optical depth ALONG THE RAY - defined as \f$ \tau=0 \f$ at the
         * front of the ray */
        double tau;
};


/** \brief Characteristic ray class.
 *
 * Holds values of \f$ s(r) \f$ and \f$ \mu(r) \f$ for each ray calculated by
 * using Bin Chen's RT formalism in Chen et al. (2007).
 *
 * */
class Characteristic
{
    public:
        /** Constructor requires a grid, otherwise our characteristics have no structure. */
        Characteristic(GridClass& grid);
        /** Save value of the geometric path length \f$ s \f$ at point \f$ k
         * \f$ along the ray. */
        void set_s(double s, unsigned int k);
        /** Save value of the direction cosine \f$ \mu \f$ at point \f$ k \f$
         * along the ray. */
        void set_mu(double mu, unsigned int k);
        /** Get the geometric coordinate \f$ s \f$ at point \f$ k \f$ along ray
         * \f$ i \f$. */
        double get_s(int k) const;
        /** Get the direction cosine \f$ \mu \f$ at point \f$ k \f$ along ray
         * \f$ i \f$. */
        double get_mu(int k) const;
        /** Return total number of points along ray. */
        unsigned int get_num_ray_pts() const;
        /** Iterator lower bound for the vector of data points along the
         * characteristic (the "back" of the ray). */
        std::vector<CharacteristicData>::const_iterator chardata_vec_begin() const;
        /** Iterator upper bound for the vector of data points along the
         * characteristic (the "front" of the ray). */
        std::vector<CharacteristicData>::const_iterator chardata_vec_end() const;

    protected:
        /** all the data along the characteristic ray */
        std::vector<CharacteristicData> m_data;
        /** Pointer to grid variables. */
        GridClass* m_grid;

    /** Lorentz factor. Needed in special relativistic flows. */
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
