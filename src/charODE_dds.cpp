#include <vector>
#include <math.h>
#include <charODE_dds.h>
#include <grid.h>

using namespace std;

charODE_dds::charODE_dds(GridClass& grid)
    : m_grid(&grid)
{}

void charODE_dds::operator() (const vector<double>& x,
                              vector<double>&       dxds,
                              const double          s)
{
    // dr/ds (Eq. 3.4a of Mihalas (1980))
    dxds.at(0) = gamma_ltz(m_grid->beta(x.at(0))) *
        (x.at(1) + m_grid->beta(x.at(0)));
    // d\mu/ds (Eq. 3.4b of Mihalas (1980))
    dxds.at(1) = gamma_ltz(m_grid->beta(x.at(0))) *
        (1.0 - pow(x.at(1), 2)) *
        (((1.0 + m_grid->beta(x.at(0)) * x.at(1)) / x.at(0)) -
         pow(gamma_ltz(m_grid->beta(x.at(0))), 2) *
         (x.at(1) + m_grid->beta(x.at(0))) *
         m_grid->dbeta_dr(x.at(0)));
}
