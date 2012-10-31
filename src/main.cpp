#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <grid.h>
#include <my_exceptions.h>
#include <characteristic.h>
#include <calc_rays.h>

using namespace std;

int main(int argc, char* argv[])
{
    cout.setf(ios::scientific);

    // user must supply YAML file as argument
    if (argc != 2) throw WrongCLIUsage();

    GridClass grid(argv[1]); // grid constructor reads YAML file

    ofstream myfile;
    myfile.open("derp.out");
    myfile.setf(ios::scientific);

    vector<Characteristic> char_ray;
    /* Initialize one fewer than the number of layers because a ray tangent to
       the outermost layer doesn't actually sample any of the atmosphere. */
    for (int i = 0; i < grid.get_num_layers()-1; i++)
    {
        Characteristic one_ray(grid, i);
        char_ray.push_back(one_ray);
    }

    // integrate characteristic ODEs
    for (vector<Characteristic>::iterator it_char = char_ray.begin();
            it_char != char_ray.end();
            it_char++)
    {
        calc_rays(grid, it_char);
    }

    // write out results
    myfile << "#mu" << "  " << "s" << "  " << "rad" << endl;
    for (vector<Characteristic>::iterator it_char = char_ray.begin();
            it_char != char_ray.end();
            it_char++)
    {
        for (int i = 0; i < 128; i++)
        {
            myfile << acos(it_char->get_mu(i)) << " "
                << it_char->get_s(i) << " "
                << grid.rad(i)
                << endl;
        }
        for (int i = 128; i < 256; i++)
        {
            myfile << acos(it_char->get_mu(i)) << " "
                << it_char->get_s(i) << " "
                << grid.rad(i-128)
                << endl;
        }
        myfile << endl;
    }

    myfile.close();

    return 0;
}
