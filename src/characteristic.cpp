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
    pair<double, double> one_pair;
    one_pair.first = s;
    one_pair.second = mu;
    m_s_mu.push_back(one_pair);
}

double Characteristic::get_s(int i)
{
    return m_s_mu.at(i).first;
}

double Characteristic::get_mu(int i)
{
    return m_s_mu.at(i).second;
}

CharNCI::CharNCI(GridClass& grid, int i)
    : Characteristic(grid),
      m_p(grid.rad(i)),
      m_tangent_layer_index(i)
{}

double CharNCI::get_p()
{
    return m_p;
}

unsigned int CharNCI::get_tangent_layer_index()
{
    return m_tangent_layer_index;
}

CharNCI_F::CharNCI_F(GridClass& grid, int i)
    : CharNCI(grid, i)
{}

double CharNCI_F::sign_of_mu()
{
    return +1.0;
}

vector< pair<double, double> >::const_iterator Characteristic::s_mu_vec_begin()
{
    return m_s_mu.begin();
}

vector< pair<double, double> >::const_iterator Characteristic::s_mu_vec_end()
{
    return m_s_mu.end();
}

void CharNCI_F::operator() (const vector<double>& s,
                            vector<double>&       dsdr,
                            const double          r)
{
    const double gamma    = gamma_ltz(m_grid->beta(r));
    const double beta     = m_grid->beta(r);
    // TODO: check if that r in the denominator of Eq. 14 in Hauschildt (1992)
    // should be there. I think it shouldn't - it keeps mu from being unitless.
    // Also it makes my ODEs explode because r is enormous (~10^14 cm).
    const double mu_E     = sign_of_mu() * sqrt(1.0 - (pow(m_p, 2) / pow(r, 2)));
    const double mu       = (mu_E - beta) / (1.0 - beta * mu_E);

    /* ds/dr, which is just the inverse of dr/ds (Eq. 5a of Hauschildt (1992)) */
    dsdr.at(0) = 1.0 / (gamma * (mu + beta));
}

CharNCI_B::CharNCI_B(GridClass& grid, int i)
    : CharNCI(grid, i)
{}

void CharNCI_B::operator() (const vector<double>& x,
                              vector<double>&       dsdr,
                              const double          r)
{
    const double gamma    = gamma_ltz(m_grid->beta(r));
    const double beta     = m_grid->beta(r);
    const double mu_E     = sign_of_mu()*sqrt(1.0 - (pow(m_p, 2) / pow(r, 2))) / r;
    const double mu       = (mu_E - beta) / (1.0 - beta * mu_E);

    /* ds/dr, which is just the inverse of dr/ds (Eq. 5a of Hauschildt (1992)) */
    dsdr.at(0) = 1.0 / (gamma * (mu + beta));
}

double CharNCI_B::sign_of_mu()
{
    return -1.0;
}

// TODO: implement () operator for odeint on core-intersecting rays
