#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_const_cgsm.h>
#include <my_exceptions.hpp>
#include <characteristic.hpp>
#include <calc_rays.hpp>
#include <misc.hpp>
#include <params.hpp>
#include <read_params.hpp>

using namespace std;

int main(int argc, char* argv[])
{
    cout << endl << "Welcome to spark!" << endl << endl;
    cout.setf(ios::scientific);

    /* User must supply YAML file as argument */
    if (argc != 2) throw WrongCLIUsage();

    /* Create data structre to hold all runtime parameters. */
    ParamsClass params;

    /* Read data in YAML file. */
    string yaml_file = string(argv[1]);
    read_params(yaml_file, &params);

    /* Now read grid data. */
    cout << "Reading in grid variables..." << endl << endl;
    GridClass grid(params.layer_file);

    /* Sanity check: velocity field must be monotonic. */
    for (unsigned int i = 0; i < grid.get_num_layers()-1; i++)
    {
        if (grid.rad(i) > grid.rad(i+1) || grid.vel(i) > grid.vel(i+1))
            throw NonmonotonicVelocityField(grid.vel(i));
    }

    // Iterate over raidus-velocity coordinate pairs. (first = radius; second = velocity)
    for (vector< pair<double, double> >::const_iterator it_grid = grid.begin(); it_grid != grid.end(); ++it_grid)
    {
        // If velocity is negative then die.
        if (it_grid->second < 0.0)
        {
            throw NegativeVelocity(it_grid->second);
        }
        // If velocity is faster than the speed of light then die.
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
    for (unsigned int i = 0; i < grid.get_num_layers()-1; i++)
    {
        CharNCI_F one_ray(grid, i);
        char_ray_front.push_back(one_ray);
    }
    for (unsigned int i = 0; i < grid.get_num_layers()-1; i++)
    {
        CharNCI_B one_ray(grid, i);
        char_ray_back.push_back(one_ray);
    }

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
