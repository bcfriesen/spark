#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <grid.h>
#include <my_exceptions.h>
#include <characteristic.h>
#include <calc_rays.h>
#include <misc.h>
#include <utils/make_linelist_hdf.h>

using namespace std;

int main(int argc, char* argv[])
{
    cout << endl << "Welcome to SPARK!" << endl << endl;
    cout.setf(ios::scientific);

    /* User must supply YAML file as argument */
    if (argc != 2) throw WrongCLIUsage();

    /* Grid constructor reads YAML file. */
    cout << "Reading in grid variables..." << endl << endl;
    GridClass grid(argv[1]);

    /* Sanity check: velocity field must be monotonic. */
    for (int i = 0; i < grid.get_num_layers()-1; i++)
    {
        if (grid.rad(i) > grid.rad(i+1) || grid.vel(i) > grid.vel(i+1))
            throw NonmonotonicVelocityField(grid.vel(i));
    }

    ofstream myfile;
    myfile.open("derp.out");
    myfile.setf(ios::scientific);

    /* Split up characteristics according to which half of the ejecta they
     * trace ("front" means "closer to observer"). By construction, rays which
     * don't intersect the core match at s=0. */
    vector<Characteristic> char_ray_front;
    vector<Characteristic> char_ray_back;

    /* Initialize non-core-intersecting rays. There are one fewer of these than
     * the number of layers because a ray tangent to the outermost layer
     * doesn't actually sample any of the atmosphere. */
    cout << "Setting up characteristic rays..." << endl << endl;
    for (int i = 0; i < grid.get_num_layers()-1; i++)
    {
        Characteristic one_ray(grid, i);
        char_ray_front.push_back(one_ray);
        char_ray_back.push_back(one_ray);
    }

    int result = make_linelist_hdf("/home5/baker/friesen/SPARK/src/utils/gf0100.all.txt");
    if (result != 0) return -1;

    /* Integrate characteristic ODEs forward from s=0. */
    cout << "Integrating forward characteristics..." << endl << endl;
    for (vector<Characteristic>::iterator it_char = char_ray_front.begin();
            it_char != char_ray_front.end();
            it_char++)
    {
        calc_rays(grid, it_char, FORWARD);
    }
    /* Integrate characteristic ODEs backward from s=0. */
    cout << "Integrating backward characteristics..." << endl << endl;
    for (vector<Characteristic>::iterator it_char = char_ray_back.begin();
            it_char != char_ray_back.end();
            it_char++)
    {
        calc_rays(grid, it_char, BACKWARD);
    }

    /* Write out results. */
    cout << "Saving characteristic data results..." << endl << endl;
    for (int i = 0; i < grid.get_num_layers()-1; i++)
    {
        myfile << grid.rad(i);
        for (vector<Characteristic>::iterator it_char = char_ray_front.begin();
                it_char != char_ray_front.end();
                it_char++)
        {
            myfile << " " << acos(it_char->get_mu(i));
        }
        myfile << endl;
    }

    myfile.close();

    return 0;
}
