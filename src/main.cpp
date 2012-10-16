#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <misc.h>
#include <grid.h>
#include <my_exceptions.h>
#include <characteristic.h>
#include <calc_rays.h>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 2) throw WrongCLIUsageException();

    GridClass grid(argv[1]);

    vector< vector<double> > x_vec;
    vector<double> times;

    cout.setf(ios::scientific);

    ofstream myfile;
    myfile.open("derp.out");

    vector<Characteristic> char_ray;
    for (int i = 0; i < grid.get_num_layers()-1; i++)
    {
        Characteristic one_ray(grid, i);
        char_ray.push_back(one_ray);
    }

    calc_rays(grid, char_ray);

    myfile.close();

    return 0;
}
