#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>

using namespace std;

//! Error message namespace.
/**
 * Shamelessly stolen from yaml-cpp.
 */
namespace ErrorMsg
{
    const string const FILE_NOT_FOUND = "File not found: ";
    const string const INTERP_OUT_OF_RANGE = "Interpolation out of range";
    const string const WRONG_CLI_USAGE = "Usage: <RTCPP_binary> <yaml_file>";
}


//! Generic exception class.
/**
 * Shamelessly stolen from yaml-cpp. Confusingly, STL contains already an
 * exception class called "exception", from which the class "runtime_error"
 * derives, from which this custom class "Exception" derives. I feel like I
 * need to draw a phylogenetic tree for the exception lineage.
 */
class Exception: public std::runtime_error
{
    public:
        Exception(const string& msg_)
            : runtime_error(build_what(msg_)),
              msg(msg_)
        {}
        virtual ~Exception() throw() {}
        std::string msg;

    private:
        static const string build_what(const string& msg)
        {
            stringstream output;
            output << "ERROR: " << msg;
            return output.str();
        }
};


//! Exception for missing files.
class FileNotFoundException: public Exception
{
    public:
        FileNotFoundException(string filename)
            : Exception(ErrorMsg::FILE_NOT_FOUND + filename)
        {}
};

//! Exception for out-of-bounds interpolation.
class InterpOutOfRangeException: public Exception
{
    public:
        InterpOutOfRangeException()
            : Exception(ErrorMsg::INTERP_OUT_OF_RANGE)
        {}
};

//! Exception for wrong # of command-line arguments.
class WrongCLIUsageException: public Exception
{
    public:
        WrongCLIUsageException()
            : Exception(ErrorMsg::WRONG_CLI_USAGE)
        {}
};
