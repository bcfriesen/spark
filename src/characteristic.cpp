#include <characteristic.h>

CharacteristicClass::CharacteristicClass(double p_in)
{
    p = p_in;
}

CharacteristicClass::~CharacteristicClass()
{

}

double CharacteristicClass::get_p()
{
    return p;
}

double CharacteristicClass::get_r(double s)
{
    // TODO: implement this
    return 0.0;
}
