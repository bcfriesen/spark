#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <vector>

using namespace std;

class GridClass
{
    public:
        GridClass();
        ~GridClass();
        void initialize(string filename);
        double get_rad(int layer) const;
        double get_vel(int layer) const;
        double get_beta(int layer) const;
    private:
        vector<double> rad;
        vector<double> vel;
        vector<double> beta;
        vector<double> dbeta_dr;
};

#endif
