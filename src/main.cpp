#include <iostream>
#include <grid.h>
#include <misc.h>

using namespace std;

int main(int argc, char* argv[])
{
    GridClass grid;
    grid.initialize("/Users/brian/RTCPP/layer.pah_std_d3.3.dat");

    // testing
    cout.setf(ios::scientific);
    cout << grid.dbeta_dr(127) << endl;
    cout << grid.vel(10) << endl;
    cout << grid.vel(11) << endl;
    cout << grid.vel((grid.rad(10) + grid.rad(11))/2.0) << endl;
    cout << grid.vel((grid.rad(10) + grid.rad(11))/1.8) << endl;

    return 0;
}
