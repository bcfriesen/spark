#include <vector>
#include <characteristic.h>

using namespace std;

Characteristic::Characteristic(GridClass& grid, int i)
    : m_grid(&grid),
      m_p(grid.rad(i)),
      tangent_layer_index(i)
{}

double Characteristic::get_p()
{
    return m_p;
}

void Characteristic::push_s(double s_)
{
    s.push_back(s_);
}

void Characteristic::push_mu(double mu_)
{
    mu.push_back(mu_);
}

double Characteristic::get_s(int i)
{
    return s.at(i);
}

double Characteristic::get_mu(int i)
{
    return mu.at(i);
}

int Characteristic::get_tangent_layer_index()
{
    return tangent_layer_index;
}