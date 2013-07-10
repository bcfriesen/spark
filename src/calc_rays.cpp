#include <iostream>
#include <cmath>

#include <grid.hpp>
#include <characteristic.hpp>

using namespace std;
void calc_rays(GridClass* grid, vector<TangentRay>& ray_vector)
{
    int i, j, k;
    i = 0;
    k = 1;
    if (k <= i+1) {
        j = k;
    } else {
        j = (2*i + 2) - k;
    }

    cout << "i = " << i << endl;
    cout << "k = " << k << endl;
    cout << "j = " << j << endl;

    const double mu_p = +sqrt(pow(grid->rad(j), 2) - pow(grid->rad(i+1), 2)) / grid->rad(j);
    const double mu_m = -sqrt(pow(grid->rad(j), 2) - pow(grid->rad(i+1), 2)) / grid->rad(j);

    cout << grid->rad(0) << endl;
    cout << grid->rad(1) << endl;
    cout << "mu_p " << mu_p << endl;
    cout << "mu_m " << mu_m << endl;

    ray_vector.at(i).set_mu(mu_p, k);
}
