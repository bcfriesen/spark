#include <iostream>
#include <cmath>
#include <vector>
#include <misc.hpp>
#include <characteristic.hpp>

using namespace std;

Characteristic::Characteristic(GridClass& grid)
    : m_grid(&grid)
{}

void Characteristic::set_s(double s, unsigned int k)
{
    m_data.at(k).s = s;
}

void Characteristic::set_mu(double mu, unsigned int k)
{
    m_data.at(k).mu = mu;
}

double Characteristic::get_s(int k) const
{
    return m_data.at(k).s;
}

double Characteristic::get_mu(int k) const
{
    return m_data.at(k).mu;
}

unsigned int Characteristic::get_num_ray_pts() const
{
    /* The vector of s points and the vector of mu points should be the same
     * size, so just pick one. */
    return m_data.size();
}

std::vector<CharacteristicData>::const_iterator Characteristic::chardata_vec_begin() const
{
    return m_data.begin();
}

std::vector<CharacteristicData>::const_iterator Characteristic::chardata_vec_end() const
{
    return m_data.end();
}

TangentRay::TangentRay(GridClass& grid, int i)
    : Characteristic(grid),
      m_p(grid.rad(i)),
      m_tangent_layer_index(i)
{
    m_data.resize(2*i + 3);
}

double TangentRay::get_p() const
{
    return m_p;
}

unsigned int TangentRay::tangent_layer_index() const
{
    return m_tangent_layer_index;
}
