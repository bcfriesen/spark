#include <iostream>
#include <cmath>
#include <vector>
#include <misc.hpp>
#include <characteristic.hpp>

using namespace std;

Characteristic::Characteristic(GridClass& grid)
    : m_grid(&grid)
{}

void Characteristic::push_s_mu(double s, double mu)
{
    m_s.push_back(s);
    m_mu.push_back(mu);
}

double Characteristic::get_s(int k) const
{
    return m_s.at(k);
}

double Characteristic::get_mu(int k) const
{
    return m_mu.at(k);
}

unsigned int Characteristic::get_num_ray_pts() const
{
    /* The vector of s points and the vector of mu points should be the same
     * size, so just pick one. */
    return m_s.size();
}

vector<double>::const_iterator Characteristic::s_vec_begin() const
{
    return m_s.begin();
}

vector<double>::const_iterator Characteristic::s_vec_end() const
{
    return m_s.end();
}

vector<double>::const_iterator Characteristic::mu_vec_begin() const
{
    return m_mu.begin();
}

vector<double>::const_iterator Characteristic::mu_vec_end() const
{
    return m_mu.end();
}

TangentRay::TangentRay(GridClass& grid, int i)
    : Characteristic(grid),
      m_p(grid.rad(i)),
      m_tangent_layer_index(i)
{}

double TangentRay::get_p() const
{
    return m_p;
}

unsigned int TangentRay::tangent_layer_index() const
{
    return m_tangent_layer_index;
}
