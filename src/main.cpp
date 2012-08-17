#include <iostream>
#include <grid.h>
#include <misc.h>

using namespace std;

int main(int argc, char* argv[])
{
    GridClass grid;
    grid.initialize("/home/friesen/RTCPP/layer.pah_std_d3.3.dat");

    // testing
    cout.setf(ios::scientific);

    return 0;
}
