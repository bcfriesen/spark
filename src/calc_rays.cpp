#include <calc_rays.h>
#include <vector>
#include <iostream>
#include <charODE_dds.h>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

void calc_rays(GridClass &grid, vector<Characteristic> &char_ray)
{

    vector<double> x(2);

    /** Upper limit of integration for the characteristic ray ODEs as functions
      * of s. We'll start this as a small fraction of the smallest radial point
      * in the grid, then increment it until r(s) reaches that radial point.
      * Then we invert the ODEs to write them as functions of r and continue
      * the integration all the way to the outermost radial point. See
      * discussion on p. 581 of Mihalas, ApJ, 237, 574 (1980). */
    double s_stop;

    /** Integrate the d/ds equations toward larger and larger s until r(s) reaches
      * r_min, the smallest radial point on the grid. Then invert. This variable
      * is the amount by which s is incremented until we reach r = r_min. */
    const double s_stop_increment = 0.01 * grid.rad(0);

    // solve the characteristic equations along each ray for s(r) and mu(r)
    for (vector<Characteristic>::iterator it_char = char_ray.begin();
         it_char != char_ray.end();
         it_char++)
    {
        cout << "working on p = " << it_char->get_p() << endl;
        s_stop = 0.0;

        // initial conditions for d/dr ODEs
        vector<double> x_ddr;

        // Calculate r(s) for gradually increasing s until r = smallest radial
        // point on grid.
        do
        {
            s_stop += s_stop_increment;

            charODE_dds ode_dds(grid);
            // TODO: implement charODE_ddr and uncomment next line
            // charODE_ddr ode_ddr;

            // initial conditions
            x.at(0) = it_char->get_p();
            x.at(1) = -grid.beta(it_char->get_p());
            // TODO: figure out a more automagic initial step size
            integrate(ode_dds, x, 0.0, s_stop, 1.0e-3*s_stop);

          // TODO: interpolate between the two values of r(s) that straddle r_min?
        } while (x.at(0) < it_char->get_p());

        // for each ray we'll store the full s(r) and mu(r) information, then
        // copy the whole vector to the appropriate ray's member vector
        vector<double> s_tmp;
        vector<double> mu_tmp;

        s_tmp.push_back(s_stop);
        mu_tmp.push_back(x.at(1));

        cout << "r = " << it_char->get_p() << endl;
        cout << "s = " << s_stop << endl;
        cout << "mu = " << x.at(1) << endl;

        // save these values of r and mu since they'll be the initial conditions
        // for each integration of the d/dr ODEs
        x_ddr = x;

        // Now that we have s(r_0) and mu(r_0), switch to the new ODEs and
        // integrate these to find s(r) and mu(r) at each radial point, using
        // s(r_0) and mu(r_0) as initial conditions.

        // TODO: implement d/dr ODEs

        // TODO: figure out how to do the whole s integration in one sweep, or
        // at least do this two-step junk under the covers
    }

}
