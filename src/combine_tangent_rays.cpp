#include <iostream>
#include <algorithm>
#include <characteristic.hpp>

using namespace std;
CharNCI combine_tangent_rays(const CharNCI_B back_ray, const CharNCI_F front_ray, GridClass& grid)
{
    if (back_ray.get_tangent_layer_index() != front_ray.get_tangent_layer_index())
    {
       // TODO: make an exception for this
       cout << "RAYS DON'T MATCH" << endl;
    }

    const unsigned int tangent_layer_index = back_ray.get_tangent_layer_index();

    CharNCI complete_ray(grid, tangent_layer_index);

    /* When we discretize a ray on a grid, it's most convenient to think about
     * the first index on the ray being at \tau = 0, i.e., where the ray
     * emerges from the atmosphere toward the observer. The index increases as
     * we trace "backwards" along the ray, through the atmosphere, until we
     * reach the side farthest from the observer.
     *
     * In this code both halves of the ray start at the center of the
     * atmosphere and move away from it. So to combine the two halves to form a
     * single ray, we first flip the data points on the front half around,
     * since those indices increase as we move *toward* the observer, and we
     * want them to increase moving *away* from the observer. The data points
     * on the back half of the ray stay unchanged, since they increase moving
     * away from the observer. */

    /* Flipping the front half over means we also need to recalculate s, the
     * geometric length along the ray. We want it to start at 0 at the front
     * edge, not at the center of the atmosphere. */
    vector<double> tmp_s(front_ray.get_num_ray_pts()), tmp_s1;
    vector<double> tmp_mu;
    for (vector< pair<double, double> >::const_iterator it_s_mu = front_ray.s_mu_vec_end(); it_s_mu != front_ray.s_mu_vec_begin(); --it_s_mu)
    {
        tmp_s1.push_back(it_s_mu->first);
        tmp_mu.push_back(it_s_mu->second);
    }
    tmp_s.at(0) = 0.0;
    for (unsigned int i = 1; i < tmp_s1.size(); ++i)
    {
        tmp_s.at(i) = tmp_s.at(i-1) + (tmp_s1.at(i-1) - tmp_s1.at(i));
    }

    for (unsigned int i = 0; i < front_ray.get_num_ray_pts(); ++i)
    {
        complete_ray.push_s_mu(tmp_s.at(i), tmp_mu.at(i));
    }

    for (unsigned int i = 0; i < complete_ray.get_num_ray_pts(); ++i)
    {
        cout << complete_ray.get_s(i) << " " << complete_ray.get_mu(i) << endl;
    }
    cout << complete_ray.get_num_ray_pts() << endl;
    cout << tmp_mu.size() << endl;

}
