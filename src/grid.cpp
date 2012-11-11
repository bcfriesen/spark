#include <fstream>
#include <algorithm>
#include <yaml-cpp/yaml.h>
#include <grid.h>
#include <my_exceptions.h>
#include <const.h>
#include <misc.h>

using namespace std;

GridClass::GridClass(char* yaml_file)
{
    /* YAML file has all the file locations, parameters, etc. */
    ifstream infile;
    infile.open(yaml_file);

    if (!infile) throw FileNotFound(yaml_file);

    string layer_file;

    try
    {
        YAML::Parser parser(infile);
        YAML::Node doc;
        parser.GetNextDocument(doc);
        // we read PHOENIX-style layer files
        doc["layer_file"] >> layer_file;
        doc["num_core_intersect_rays"] >> num_core_intersect_rays;
    }
    catch(YAML::ParserException& e)
    {
        cout << e.what() << endl;
    }

    infile.close();

    infile.open(layer_file.c_str());

    if (!infile) throw FileNotFound(layer_file);

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
}

int GridClass::get_num_layers() const
{
    return rad_vel.size();
}

int GridClass::get_num_core_intersect_rays() const
{
    return num_core_intersect_rays;
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
    return rad_vel.at(layer).second / physconst::c_light;
}

double GridClass::beta(double rad) const
{
    return vel(rad) / physconst::c_light;
}

double GridClass::dbeta_dr(int layer) const
{
    /* Forward derivative at the back. */
    if (layer == 0)
    {
        return (beta(layer+1) - beta(layer)) /
            (rad_vel.at(layer+1).first - rad_vel.at(layer).first);
    }
    /* Backward derivative at the front. */
    else if (layer == rad_vel.size()-1)
    {
        return (beta(layer) - beta(layer-1)) /
            (rad_vel.at(layer).first - rad_vel.at(layer-1).first);
    }
    /* Centered derivative anywhere in between. */
    else
    {
        return (beta(layer+1) - beta(layer-1)) /
            (rad_vel.at(layer+1).first - rad_vel.at(layer-1).first);
    }
}

double GridClass::dbeta_dr(double rad) const
{
    // TODO: figure out how to make step size not so knobby
    const double dr = 1.0e-4 * rad_vel.at(0).first;
    /* Forward derivative at the back. */
    if (rad < rad_vel.at(1).first)
    {
        return (beta(rad+dr) - beta(rad)) / dr;
    }
    /* Backward derivative at the front. */
    else if (rad > rad_vel.at(rad_vel.size()-2).first)
    {
        return (beta(rad) - beta(rad-dr)) / dr;
    }
    /* Centered derivative anywhere in between. */
    else
    {
        return (beta(rad+dr) - beta(rad-dr)) / (2.0*dr);
    }
}

vector< pair<double, double> >::iterator GridClass::begin()
{
    return rad_vel.begin();
}

vector< pair<double, double> >::iterator GridClass::end()
{
    return rad_vel.end();
}
