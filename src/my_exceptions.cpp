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

FileNotFoundException::FileNotFoundException(string filename)
    : Exception(ErrorMsg::FILE_NOT_FOUND + filename)
{}

InterpOutOfRangeException::InterpOutOfRangeException(double x)
: Exception(ErrorMsg::INTERP_OUT_OF_RANGE + boost::lexical_cast<string>(x))
{}

WrongCLIUsageException::WrongCLIUsageException()
: Exception(ErrorMsg::WRONG_CLI_USAGE)
{}
