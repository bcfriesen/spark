#include <fstream>
#include <sstream>
#include <grid.h>
#include <griddata.h>

using namespace std;

GridClass::GridClass(string filename)
{
    ifstream infile;
    infile.open(filename.c_str());
    if (!infile)
    {
        cout << "ERROR: could not open file: " << filename << endl;
        exit(1);
    }
    cout << "Opening file: " << filename << endl;

    string oneline;
    double x1, x2, x3, x4, x5;
    GridDataClass* gd_oneline = new GridDataClass;
    getline(infile, oneline); // first line has # of layers
    istringstream is(oneline);
    int nlayer;
    is >> x1 >> nlayer;
    cout << "# of layers: " << nlayer << endl;

    for (int i = 0; i < nlayer; i++)
    {
        getline(infile, oneline);
        istringstream is(oneline);
        is >> x1 >> x2 >> x3 >> x4 >> x5;
        gd_oneline->set_rad(x1);
        gd_oneline->set_vel(x5);
        griddata.push_back(*gd_oneline);
    }

    // layers are reversed in layer file

    infile.close();
    delete gd_oneline;
}

GridClass::~GridClass() { }

double GridClass::get_rad(int layer)
{
    GridDataClass* tmp = new GridDataClass;
    *tmp = griddata.at(layer);
    double rad = tmp->get_rad();
    delete tmp;
    return rad;
}

double GridClass::get_vel(int layer)
{
    GridDataClass* tmp = new GridDataClass;
    *tmp = griddata.at(layer);
    double vel = tmp->get_vel();
    delete tmp;
    return vel;
}
