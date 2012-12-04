#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_const_cgsm.h>
#include <grid.h>
#include <my_exceptions.h>
#include <characteristic.h>
#include <calc_rays.h>
#include <misc.h>
#include <utils/make_linelist_hdf.h>

using namespace std;

int main(int argc, char* argv[])
{
    cout << endl << "Welcome to spark!" << endl << endl;
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

    // TODO: include other sanity checks (v > 0, v < c, etc.)
    for (vector< pair<double, double> >::const_iterator it_grid = grid.begin(); it_grid != grid.end(); ++it_grid)
    {
        if (it_grid->second < 0.0)
        {
            throw NegativeVelocity(it_grid->second);
        }
        else if (it_grid->second > GSL_CONST_CGSM_SPEED_OF_LIGHT)
        {
            throw SuperluminalVelocity(it_grid->second);
        }
    }

    ofstream myfile;
    // TODO: make output file name an option in the YAML file
    myfile.open("derp.out");
    myfile.setf(ios::scientific);

    /* Split up characteristics according to which half of the ejecta they
     * trace ("front" means "closer to observer"). By construction, rays which
     * don't intersect the core match at s=0. */
    vector<CharNCI_F> char_ray_front;
    vector<CharNCI_B> char_ray_back;

    /* Initialize non-core-intersecting rays. There are one fewer of these than
     * the number of layers because a ray tangent to the outermost layer
     * doesn't actually sample any of the atmosphere. */
    cout << "Setting up characteristic rays..." << endl << endl;
    for (int i = 0; i < grid.get_num_layers()-1; i++)
    {
        CharNCI_F one_ray(grid, i);
        char_ray_front.push_back(one_ray);
    }
    for (int i = 0; i < grid.get_num_layers()-1; i++)
    {
        CharNCI_B one_ray(grid, i);
        char_ray_back.push_back(one_ray);
    }

    // TODO: make linelist construction option a parameter in the YAML file
    // make_linelist_hdf("../src/utils/gfall.dat");

    /* Integrate characteristic ODEs forward from s=0. */
    cout << "Integrating forward characteristics..." << endl << endl;
    calc_rays(grid, char_ray_front);

    /* Write out results. */
    cout << "Saving characteristic data results..." << endl << endl;
    for (vector<CharNCI_F>::iterator it_char = char_ray_front.begin(); it_char != char_ray_front.end(); ++it_char)
    {
        for (vector< pair<double, double> >::const_iterator it_val = it_char->s_mu_vec_begin(); it_val != it_char->s_mu_vec_end(); ++it_val)
        {
            myfile << " " << it_val->first << " " << it_val->second << endl;
        }
        myfile << endl;
    }

    myfile.close();

    return 0;
}
