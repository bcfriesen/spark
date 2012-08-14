#ifndef CHARACTERISTIC_H
#define CHARACTERISTIC_H

class CharacteristicClass
{
    public:
        CharacteristicClass(double p_in);
        ~CharacteristicClass();
        double get_p();
        double get_r(double s);
        double get_mu(double s);
    private:
        double p; // impact parameter
};

#endif
