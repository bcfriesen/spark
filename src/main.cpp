#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <grid.h>
#include <misc.h>
#include <my_exceptions.h>
#include <characteristic.h>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

int main(int argc, char* argv[])
{
    if (argc != 2) throw WrongCLIUsageException();

    // read the YAML file and get layer data for grid
    GridClass grid(argv[1]);
    // the characteristic ODEs are unique for each ray
    vector<CharacteristicODEClass> charODE;

    // skip the outermost shell - a ray tangent to that shell never samples
    // the ejecta
    for (int i = 0; i < grid.get_num_layers()-1; i++)
    {
        // FIXME: is this many initializations horrible for I/O?
        CharacteristicODEClass one_charODE(grid, i);
        charODE.push_back(one_charODE);
    }

    vector<double> x(2);
    vector< vector<double> > x_vec;
    vector<double> times;

    cout.setf(ios::scientific);

    ofstream myfile;
    myfile.open("derp.out");

    for (vector<CharacteristicODEClass>::iterator it_r = charODE.begin(); it_r != charODE.end(); it_r++)
    {
        cout << "working on p = " << it_r->get_p() << endl;
        // set up s > 0 grid
        // TODO: figure out a better way to distribute s points
        vector<double> spts(64);
        for (int i = 0; i < spts.size(); i++)
        {
            spts.at(i) = double(i) * (1.0e15 / double(spts.size()));
        }

        // calculate r(s) at each point s > 0
        for (vector<double>::iterator it_s = spts.begin(); it_s != spts.end(); it_s++)
        {
            // ODE initial conditions. see Mihalas (1980) or Hauschildt (1993)
            x.at(0) = it_r->get_p();
            x.at(1) = -grid.beta(it_r->get_p());
            try
            {
                integrate(*it_r, x, 0.0, *it_s, 1.0e-5*grid.rad(0));
                x_vec.push_back(x);
                times.push_back(*it_s);
            }
            // if we integrate past the edge of the grid, just stop integrating and
            // use whatever answer had from the last valid integration point
            catch (InterpOutOfRangeException& e) { break; }
        }
        // repeat r(s) calculation for s < 0
        // TODO: figure out how to do the whole s integration in one sweep, or
        // at least do this two-step junk under the covers
        // TODO: figure out a better way to distribute s points
        for (int i = 0; i < spts.size(); i++)
        {
            spts.at(i) = double(i) * (-1.0e15 / double(spts.size()));
        }
        // calculate r(s) at each point
        for (vector<double>::iterator it_s = spts.begin(); it_s != spts.end(); it_s++)
        {
            x.at(0) = it_r->get_p();
            x.at(1) = -grid.beta(it_r->get_p());
            try
            {
                integrate(*it_r, x, 0.0, *it_s, -1.0e-5*grid.rad(0));
                x_vec.push_back(x);
                times.push_back(*it_s);
            }
            catch (InterpOutOfRangeException& e) { break; }
        }

        for (int i = 0; i < x_vec.size(); i++)
        {
            myfile << times.at(i) << "\t"
                   << x_vec.at(i).at(0) << "\t"
                   << x_vec.at(i).at(1) << "\t"
                   << acos(x_vec.at(i).at(1)) << "\t"
                   << endl;
        }
    }

    myfile.close();


    return 0;
}
