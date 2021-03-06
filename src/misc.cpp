#include <limits>
#include <cmath>
#include <misc.hpp>
#include <my_exceptions.hpp>

using namespace std;

double gamma_ltz(double beta)
{
    return 1.0 / sqrt(1.0 - pow(beta, 2));
}

double interpolate(vector< pair<double, double> > table, double x)
{
    const double inf = numeric_limits<double>::infinity();

    // make sure we're interpolating within limits
    if (x - table.back().first > 0.0 || x - table.at(0).first < 0.0)
        throw InterpOutOfRange(x);

    vector< pair< double, double > >::iterator it, it2;

    it = lower_bound(table.begin(), table.end(), make_pair(x, -inf));
    if (it == table.begin()) return it->second;
    it2 = it;
    --it2;
    return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}

std::string getEnvVar(std::string const &key)
{
    char* val = getenv(key.c_str());
    return val == NULL ? std::string("") : std::string(val);
}

double gaussian(double x, double a, double b, double c)
{
    return (a * exp(-pow(x - b, 2) / (2.0 * pow(c, 2))));
}
