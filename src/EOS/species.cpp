#include <element.h>

void SpeciesClass::set_abund(double N)
{
    m_N = N;
}

double SpeciesClass::get_abund() const
{
    return m_N;
}
