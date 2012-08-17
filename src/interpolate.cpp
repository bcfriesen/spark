// stole this whole routine from Daniel Fleischman on stackoverflow

#include <limits>
#include <vector>
#include <algorithm>
#include <iostream>
#include <my_exceptions.h>
#include <interpolate.h>

using namespace std;

double interpolate(vector<pair <double, double> > table, double x)
{
    const double inf = numeric_limits<double>::infinity();

    // make sure we're interpolating within limits
    try
    {
        double check = x - table.back().first;
        if (check > 0) throw InterpOutOfRangeException(table, x);
    }
    catch (InterpOutOfRangeException& ior)
    {
        cout << ior.what() << endl;
        // don't need to return a value but I don't know what else to do here
        return -inf;
    }
    try
    {
        double check = x - table.at(0).first;
        if (check < 0) throw InterpOutOfRangeException(table, x);
    }
    catch (InterpOutOfRangeException& ior)
    {
        cout << ior.what() << endl;
        return -inf;
    }

    vector< pair< double, double > >::iterator it, it2;

    it = lower_bound(table.begin(), table.end(), make_pair(x, -inf));
    if (it == table.begin()) return it->second;
    it2 = it;
    --it2;
    return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}
