#include <fstream>
#include <sstream>
#include <algorithm>
#include <grid.h>
#include <misc.h>
#include <my_exceptions.h>

using namespace std;

GridClass::GridClass() { }

GridClass::~GridClass() { }

void GridClass::initialize(string filename)
{
    const double c_light = 2.99792458e10;
    const double cm2km = 1.0e-5;
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

    for (int i = 0; i < nlayer; i++)
    {
        getline(infile, oneline);
        istringstream is(oneline);
        is >> x1 >> x2 >> x3 >> x4 >> x5;
        rad.push_back(x1);
        vel.push_back(x5/cm2km);
    }

    infile.close();

    // layers are reversed in the layer file, so reverse them back
    reverse(rad.begin(), rad.end());
    reverse(vel.begin(), vel.end());

    // calculate Lorentz factor beta (= v/c) for each layer
    for (vector<double>::iterator it = vel.begin(); it != vel.end(); it++)
    {
        beta.push_back(*it / c_light);
    }

    dbeta_dr.resize(beta.size());
    for (int i = 0; i < beta.size(); i++)
    {
        if (i == 0)
        {
            dbeta_dr.at(i) = (beta.at(i+1) - beta.at(i)) / (rad.at(i+1) - rad.at(i));
        }
        else if (i == beta.size()-1)
        {
            dbeta_dr.at(i) = (beta.at(i) - beta.at(i-1)) / (rad.at(i) - rad.at(i-1));
        }
        else
        {
            dbeta_dr.at(i) = (beta.at(i+1) - beta.at(i-1)) / (rad.at(i+1) - rad.at(i-1));
        }
    }

}

double GridClass::get_rad(int layer) const
{
    return rad.at(layer);
}

double GridClass::get_vel(int layer) const
{
    return vel.at(layer);
}

double GridClass::get_beta(int layer) const
{
    return beta.at(layer);
}
