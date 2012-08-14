#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <vector>
#include <griddata.h>

using namespace std;

class GridClass
{
    public:
        GridClass(string filename);
        ~GridClass();
        double get_rad(int layer);
        double get_vel(int layer);
    private:
        vector<GridDataClass> griddata;
};

#endif
