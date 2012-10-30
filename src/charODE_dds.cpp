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
    const double r        = x.at(0);
    const double mu       = x.at(1);
    const double gamma    = gamma_ltz(m_grid->beta(r));
    const double beta     = m_grid->beta(r);
    const double dbeta_dr = m_grid->dbeta_dr(r);

    double dr_ds;
    double dmu_ds;

    // dr/ds (Eq. 3.4a of Mihalas (1980))
    dr_ds = gamma * (mu + beta);

    // d\mu/ds (Eq. 3.4b of Mihalas (1980))
    dmu_ds = gamma * (1.0 - pow(mu, 2)) *
        (((1.0 + beta * mu) / r) -
         pow(gamma, 2) * (mu + beta) * dbeta_dr);

    dxds.at(0) = dr_ds;
    dxds.at(1) = dmu_ds;
}
