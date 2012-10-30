#include <vector>
#include <math.h>
#include <characteristic.h>

using namespace std;

Characteristic::Characteristic(GridClass& grid, int i)
    : m_grid(&grid),
      m_p(grid.rad(i))
{}

//! get impact parameter for this ray (at v = 0)
double Characteristic::get_p()
{
    return m_p;
}
