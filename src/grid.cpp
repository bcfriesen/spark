#include <fstream>
#include <sstream>
#include <algorithm>
#include <grid.h>
#include <misc.h>
#include <my_exceptions.h>
#include <const.h>
#include <interpolate.h>

using namespace std;

GridClass::GridClass() { }

GridClass::~GridClass() { }

void GridClass::initialize(string filename)
{
    ifstream infile;

    // make sure file exists
    try
    {
        infile.open(filename.c_str());
        if (!infile)
        {
            throw FileNotFoundException(filename);
        }
    }
    catch (FileNotFoundException& fnf)
    {
        cout << fnf.what() << endl;
        return;
    }

    cout << "Opening layer file: " << filename << endl;

    string oneline;
    double x1, x2, x3, x4, x5;
    getline(infile, oneline); // first line has # of layers
    istringstream is(oneline);
    int nlayer;
    is >> x1 >> nlayer; // don't need first number
    cout << "# of layers: " << nlayer << endl;

    pair<double, double> rv_one; // radius/velocity of one layer

    for (int i = 0; i < nlayer; i++)
    {
        getline(infile, oneline);
        istringstream is(oneline);
        is >> x1 >> x2 >> x3 >> x4 >> x5;
        rv_one.first = x1;
        rv_one.second = x5/cm2km;
        rad_vel.push_back(rv_one);
    }

    infile.close();

    // layers are reversed in the layer file, so reverse them back
    reverse(rad_vel.begin(), rad_vel.end());
}

double GridClass::rad(int layer) const
{
    return rad_vel.at(layer).first;
}

double GridClass::vel(int layer) const
{
    return rad_vel.at(layer).second;
}

double GridClass::vel(double rad) const
{
    return interpolate(rad_vel, rad);
}

double GridClass::beta(int layer) const
{
    return rad_vel.at(layer).second / cm2km / c_light;
}

double GridClass::dbeta_dr(int layer) const
{
    if (layer == 0)
    {
        // forward derivative at the back
        return (beta(layer+1) - beta(layer)) / (rad_vel.at(layer+1).first - rad_vel.at(layer).first);
    }
    else if (layer == rad_vel.size()-1)
    {
        // backward derivative at the front
        return (beta(layer) - beta(layer-1)) / (rad_vel.at(layer).first - rad_vel.at(layer-1).first);
    }
    else
    {
        // centered derivative anywhere in between
        return (beta(layer+1) - beta(layer-1)) / (rad_vel.at(layer+1).first - rad_vel.at(layer-1).first);
    }
}

