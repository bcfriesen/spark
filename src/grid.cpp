#include <fstream>
#include <sstream>
#include <algorithm>
#include <grid.h>
#include <misc.h>
#include <my_exceptions.h>
#include <const.h>
#include <interpolate.h>
#include <yaml-cpp/yaml.h>

using namespace std;

GridClass::GridClass(char* yaml_file)
{
    ifstream infile;
    infile.open(yaml_file);
    string layer_file;

    try
    {
        YAML::Parser parser(infile);
        YAML::Node doc;
        parser.GetNextDocument(doc);
        doc["layer_file"] >> layer_file;
    }
    catch(YAML::ParserException& e)
    {
        cout << e.what() << endl;
    }

    infile.close();

    infile.open(layer_file.c_str());
    if (!infile) throw FileNotFoundException(layer_file);

    string oneline;
    double x1, x2, x3, x4, x5;
    getline(infile, oneline); // first line has # of layers
    istringstream is(oneline);
    int nlayer;
    is >> x1 >> nlayer; // don't need first number

    pair<double, double> rv_one; // radius & velocity of one layer

    for (int i = 0; i < nlayer; i++)
    {
        getline(infile, oneline);
        istringstream is(oneline);
        is >> x1 >> x2 >> x3 >> x4 >> x5;
        rv_one.first = x1;
        rv_one.second = x5 / physconst::cm2km;
        rad_vel.push_back(rv_one);
    }

    infile.close();

    // sort layers with increasing radius
    sort(rad_vel.begin(), rad_vel.end());
}

GridClass::~GridClass() { }

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
    return rad_vel.at(layer).second / physconst::cm2km / physconst::c_light;
}

double GridClass::dbeta_dr(int layer) const
{
    if (layer == 0)
    {
        // forward derivative at the back
        return (beta(layer+1) - beta(layer)) /
            (rad_vel.at(layer+1).first - rad_vel.at(layer).first);
    }
    else if (layer == rad_vel.size()-1)
    {
        // backward derivative at the front
        return (beta(layer) - beta(layer-1)) /
            (rad_vel.at(layer).first - rad_vel.at(layer-1).first);
    }
    else
    {
        // centered derivative anywhere in between
        return (beta(layer+1) - beta(layer-1)) /
            (rad_vel.at(layer+1).first - rad_vel.at(layer-1).first);
    }
}

