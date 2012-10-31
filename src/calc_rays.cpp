#include <vector>
#include <iostream>
#include <cmath>
#include <calc_rays.h>
#include <charODE_dds.h>
#include <charODE_dds_new.h>
#include <charODE_ddr.h>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

void calc_rays(GridClass &grid, std::vector<Characteristic>::iterator it_char)
{
    // r and mu coordinates for d/ds equations
    vector<double> x_dds(2);
    // s and mu coordinates for d/dr equations
    vector<double> x_ddr(2);

    /** Upper limit of integration for the d/ds characteristic ray ODEs. We
     * integrate forward a tiny bit in s, then invert the ODEs to write them as
     * functions of r and continue the integration all the way to the outermost
     * radial point. (See discussion on p. 581 of Mihalas, ApJ, 237, 574 (1980)).
     * There's nothing special about the value I picked here: I just multiplied
     * it by some fraction of the smallest radial grid point so that 1.) it will
     * be on the same scale as the rest of the grid so as to avoid machine
     * roundoff errors, and 2.) it will still be cheap to integrate because it's
     * not far from the initial conditions. */
    const double s_stop = 1.0e-3 * grid.rad(0);

    for (int i = 0; i < grid.get_num_layers(); i++)
    {
        const double r = grid.rad(i);
        const double p = it_char->get_p();
        const double beta = grid.beta(i);
        // direction cosines in Eulerian frame, toward the observer (m) and away
        // from the observer (p)
        const double mu_E_p = sqrt(1.0 - pow(p / r, 2)) / r;
        const double mu_E_m = -mu_E_p;
        // direction cosines in Lagrangian frame
        const double mu_L_p = (mu_E_p - beta) / (1.0 - beta * mu_E_p);
        const double mu_L_m = (mu_E_m - beta) / (1.0 - beta * mu_E_m);
        // set up d/ds ODEs
        charODE_dds_new ode_dds_new(grid, mu_L_p);
        double r_of_s = p;
        do
        {
            integrate(ode_dds_new, r_of_s, 1.0e-6*s_stop, s_stop, 1.0e-5*s_stop);
            cout << "r is now " << r_of_s << endl;
        } while (r_of_s <= r);
    }

    //---------------- FIRST PASS: POSITIVE S -----------------
    // initial conditions for r and mu (see Eq. 3.6 and the paragraph preceding
    // it in Mihalas (1980))
    x_dds.at(0) = it_char->get_p();
    x_dds.at(1) = -grid.beta(it_char->get_p());

    // set up d/ds ODEs
    charODE_dds ode_dds(grid);

    // Now integrate d/ds equations. First do it toward positive s, then
    // again toward negative s.
    // TODO: figure out a more automagic initial step size
    integrate(ode_dds, x_dds, 1.0e-6*s_stop, s_stop, 1.0e-5*s_stop);

    // The results from this integration (r(s) and mu(s)) are now the
    // initial conditions for the d/dr equations.
    charODE_ddr ode_ddr(grid);

    // initial conditions for d/dr ODEs
    x_ddr.at(0) = s_stop;
    x_ddr.at(1) = x_dds.at(1);
    cout << "initial conditions for d/dr sweep after positive s: ";
    cout << "r small = " << x_dds.at(0) << " ";
    cout << "s(r = small) = " << x_ddr.at(0) << " ";
    cout << "mu(r = small) = " << x_ddr.at(1) << endl;

    // Now that we have s(r_0) and mu(r_0), switch to the new ODEs and
    // integrate these to find s(r) and mu(r) at each radial point, using
    // s(r_0) and mu(r_0) as initial conditions.
    integrate(ode_ddr, x_ddr, x_dds.at(0), grid.rad(0), 1.0e-5*s_stop);
    it_char->push_s(x_ddr.at(0));
    it_char->push_mu(x_ddr.at(1));
    // The stupid bookkeeping part is done. We have s(r) and mu(r) at the
    // smallest radial point on the grid. Now we just integrate from one
    // radial point to the next until we make it all the way out.
    for (int i = 0; i < grid.get_num_layers()-1; i++)
    {
        // TODO: figure out more automagic initial step size
        integrate(ode_ddr, x_ddr, grid.rad(i), grid.rad(i+1), 1.0e-5*s_stop);
        it_char->push_s(x_ddr.at(0));
        it_char->push_mu(x_ddr.at(1));
    }

    //---------------- SECOND PASS: NEGATIVE S -----------------
    // same initial conditions as before
    x_dds.at(0) = it_char->get_p();
    x_dds.at(1) = -grid.beta(it_char->get_p());
    // integrate toward negative s this time
    integrate(ode_dds, x_dds, -1.0e-6*s_stop, -s_stop, -1.0e-5*s_stop);
    // initial conditions for d/dr ODEs
    x_ddr.at(0) = -s_stop;
    x_ddr.at(1) = x_dds.at(1);
    cout << "initial conditions for d/dr sweep after negative s: ";
    cout << "r small = " << x_dds.at(0) << " ";
    cout << "s(r = small) = " << x_ddr.at(0) << " ";
    cout << "mu(r = small) = " << x_ddr.at(1) << endl;
    // integrate over same r values
    integrate(ode_ddr, x_ddr, x_dds.at(0), grid.rad(0), 1.0e-5*s_stop);
    it_char->push_s(x_ddr.at(0));
    it_char->push_mu(x_ddr.at(1));
    for (int i = 0; i < grid.get_num_layers()-1; i++)
    {
        // TODO: figure out more automagic initial step size
        integrate(ode_ddr, x_ddr, grid.rad(i), grid.rad(i+1), 1.0e-5*s_stop);
        it_char->push_s(x_ddr.at(0));
        it_char->push_mu(x_ddr.at(1));
    }

}
