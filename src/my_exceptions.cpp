#include <boost/lexical_cast.hpp>
#include <my_exceptions.h>

using namespace std;

Exception::Exception(const string& msg_)
: runtime_error(build_what(msg_)),
  msg(msg_)
{}

Exception::~Exception() throw()
{}

const string Exception::build_what(const string& msg)
{
    stringstream output;
    output << "ERROR: " << msg;
    return output.str();
}

FileNotFound::FileNotFound(string filename)
    : Exception(ErrorMsg::FILE_NOT_FOUND + filename)
{}

InterpOutOfRange::InterpOutOfRange(double x)
: Exception(ErrorMsg::INTERP_OUT_OF_RANGE + boost::lexical_cast<string>(x))
{}

WrongCLIUsage::WrongCLIUsage()
: Exception(ErrorMsg::WRONG_CLI_USAGE)
{}

NonmonotonicVelocityField::NonmonotonicVelocityField(double vel)
: Exception(ErrorMsg::NONMONOTONIC + boost::lexical_cast<string>(vel))
{}

InvalidElement::InvalidElement(int Z, int I)
: Exception(ErrorMsg::INVALID_ELEMENT + "(" + boost::lexical_cast<string>(Z) + ", " + boost::lexical_cast<string>(I) + ")")
{}
