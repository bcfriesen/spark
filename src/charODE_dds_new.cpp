#include <charODE_dds_new.h>
#include <grid.h>

using namespace std;

charODE_dds_new::charODE_dds_new(GridClass& grid, double mu)
    : m_grid(&grid),
      m_mu(mu)
{}

void charODE_dds_new::operator() (const double& r,
                                  double&       dr_ds,
                                  const double  s)
{
    const double gamma    = gamma_ltz(m_grid->beta(r));
    const double beta     = m_grid->beta(r);

    // dr/ds (Eq. 3.4a of Mihalas (1980))
    dr_ds = gamma * (m_mu + beta);
}
