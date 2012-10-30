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
    cout.setf(ios::scientific);

    // user must supply YAML file as argument
    if (argc != 2) throw WrongCLIUsageException();

    GridClass grid(argv[1]); // grid constructor reads YAML file

    ofstream myfile;
    myfile.open("derp.out");

    vector<Characteristic> char_ray;
    /* Initialize one fewer than the number of layers because a ray tangent to
       the outermost layer doesn't actually sample any of the atmosphere. */
    for (int i = 0; i < grid.get_num_layers()-1; i++)
    {
        Characteristic one_ray(grid, i);
        char_ray.push_back(one_ray);
    }

    // integrate characteristic ODEs
    calc_rays(grid, char_ray);

    myfile.close();

    return 0;
}
