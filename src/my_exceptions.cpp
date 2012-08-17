#include <my_exceptions.h>
#include <boost/lexical_cast.hpp>

FileNotFoundException::FileNotFoundException(string filename)
{
    errmsg = "File not found: " + filename;
}

FileNotFoundException::~FileNotFoundException() throw() { }

const char* FileNotFoundException::what() const throw()
{
    return this->errmsg.c_str();
}


InterpOutOfRangeException::InterpOutOfRangeException(vector< pair<double, double> > table, double x)
{
    errmsg = "Interpolation value out of range!\n";
    errmsg += "minimum: " + boost::lexical_cast<string>(table.at(0).first) + "\n";
    errmsg += "maximum: " + boost::lexical_cast<string>(table.back().first) + "\n";
    errmsg += "requested: " + boost::lexical_cast<string>(x);
}

InterpOutOfRangeException::~InterpOutOfRangeException() throw() { }

const char* InterpOutOfRangeException::what() const throw()
{
    return this->errmsg.c_str();
}
