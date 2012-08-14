#ifndef GRIDDATA_H
#define GRIDDATA_H

class GridDataClass
{
    public:
        GridDataClass();
        ~GridDataClass();
        void set_vel(double gridvel);
        void set_rad(double gridrad);
        double get_vel();
        double get_rad();
    private:
        double vel;
        double rad;
};

#endif
