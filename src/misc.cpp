#include <limits>
#include <math.h>
#include <misc.h>
#include <my_exceptions.h>

using namespace std;

double gamma_ltz(double beta)
{
    return 1.0 / sqrt(1.0 - pow(beta, 2));
}

double interpolate(vector<pair <double, double> > table, double x)
{
    const double inf = numeric_limits<double>::infinity();

    // make sure we're interpolating within limits
    if (x - table.back().first > 0 || x - table.at(0).first < 0)
        throw InterpOutOfRange(x);

    vector< pair< double, double > >::iterator it, it2;

    it = lower_bound(table.begin(), table.end(), make_pair(x, -inf));
    if (it == table.begin()) return it->second;
    it2 = it;
    --it2;
    return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}
