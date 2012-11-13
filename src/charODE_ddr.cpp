#include <vector>
#include <math.h>
#include <charODE_ddr.h>
#include <grid.h>

using namespace std;

charODE_ddr::charODE_ddr(GridClass& grid)
    : m_grid(&grid)
{}

void charODE_ddr::operator() (const vector<double>& x,
                              vector<double>&       dxdr,
                              const double          r)
{
    const double s        = x.at(0);
    const double mu       = x.at(1);
    const double gamma    = gamma_ltz(m_grid->beta(r));
    const double beta     = m_grid->beta(r);
    const double dbeta_dr = m_grid->dbeta_dr(r);

    /* ds/dr, which is just the inverse of dr/ds (Eq. 3.4a of Mihalas (1980)) */
    const double ds_dr = 1.0 / (gamma * (mu + beta));

    /* d\mu/dr, which we get from the chain rule on d\mu/ds (Eq. 3.4b of Mihalas
     * (1980)) */
    const double dmu_dr = ((1.0 - pow(mu, 2)) / (mu + beta)) *
        (((1.0 + beta * mu) / r) - pow(gamma, 2) *
         (mu + beta) * dbeta_dr);

    dxdr.at(0) = ds_dr;
    dxdr.at(1) = dmu_dr;
}
