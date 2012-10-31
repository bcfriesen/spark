#include <charODE_dds_new.h>
#include <grid.h>

using namespace std;

charODE_dds_new::charODE_dds_new(GridClass& grid, const double mu)
    : m_grid(&grid),
      m_mu(mu)
{}

void charODE_dds_new::operator() (const vector<double>& r,
                                  vector<double>&       dr_ds,
                                  const double          s)
{
    const double gamma = gamma_ltz(m_grid->beta(r.at(0)));
    const double beta  = m_grid->beta(r.at(0));

    // dr/ds (Eq. 3.4a of Mihalas (1980))
    dr_ds.at(0) = gamma * (m_mu + beta);
}
