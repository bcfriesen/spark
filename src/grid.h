#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <my_exceptions.h>
#include <const.h>
#include <interpolate.h>

using namespace std;

/**
 * Contains most grid-based variables and functions. Which is pretty much
 * everything.
 */
class GridClass
{
    public:
        //! Does most initialization.
        /**
         * Reads the radius/velocity grid data from the layer file and stores it
         * in member variables.
         */
        GridClass(char* yaml_file)
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

        //! Get number of layers in layer file.
        int get_num_layers() const { return rad_vel.size(); }

        //! Get radius from tabulated data in layer file.
        double rad(int layer) const { return rad_vel.at(layer).first; }

        //! Get velocity from tabulated data in layer file.
        double vel(int layer) const { return rad_vel.at(layer).second; }

        //! Interpolate velocity between tabulated data points.
        double vel(double rad) const { return interpolate(rad_vel, rad); }

        //!  Get \f$ v/c \f$ from tabulated data in layer file.
        double beta(int layer) const
        {
            return rad_vel.at(layer).second / physconst::c_light;
        }

        //! Interpolate \f$ v/c \f$ from tabulated data in layer file.
        double beta(double rad) const { return vel(rad) / physconst::c_light; }

        //! Calculate \f$ d\beta/dr \f$ at tabulated data points.
        double dbeta_dr(int layer) const
        {
            // forward derivative at the back
            if (layer == 0)
            {
                return (beta(layer+1) - beta(layer)) /
                    (rad_vel.at(layer+1).first - rad_vel.at(layer).first);
            }
            // backward derivative at the front
            else if (layer == rad_vel.size()-1)
            {
                return (beta(layer) - beta(layer-1)) /
                    (rad_vel.at(layer).first - rad_vel.at(layer-1).first);
            }
            // centered derivative anywhere in between
            else
            {
                return (beta(layer+1) - beta(layer-1)) /
                    (rad_vel.at(layer+1).first - rad_vel.at(layer-1).first);
            }
        }

        //! Calculate \f$ d\beta/dr \f$ by interpolating between data points.
        // FIXME: this whole thing is probably completely wrong
        double dbeta_dr(double rad) const
        {
            // TODO: figure out how to make step size not so knobby
            const double dr = 1.0e-4 * rad_vel.at(0).first;
            // forward derivative at the back
            if (rad < rad_vel.at(1).first)
            {
                return (beta(rad+dr) - beta(rad)) / dr;
            }
            // backward derivative at the front
            else if (rad > rad_vel.at(rad_vel.size()-2).first)
            {
                return (beta(rad) - beta(rad-dr)) / dr;
            }
            // centered derivative anywhere in between
            else
            {
                return (beta(rad+dr) - beta(rad-dr)) / (2.0*dr);
            }
        }

        friend double interpolate(vector< pair<double, double> > table, double x);
        friend double gamma_ltz(double beta);

    private:
        //! radius and velocity of each layer
        vector< pair<double, double> > rad_vel;
};

#endif
