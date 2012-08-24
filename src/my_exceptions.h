#ifndef MY_EXCEPTIONS_H
#define MY_EXCEPTIONS_H

#include <string>
#include <stdexcept>

//! Error message namespace.
/**
 * Shamelessly stolen from yaml-cpp.
 */
namespace ErrorMsg
{
    const std::string const FILE_NOT_FOUND = "File not found: ";
    const std::string const INTERP_OUT_OF_RANGE = "Interpolation out of range for value: ";
    const std::string const WRONG_CLI_USAGE = "Usage: <RTCPP_binary> <yaml_file>";
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
        Exception(const std::string& msg_);
        virtual ~Exception() throw();
        std::string msg;
    private:
        static const std::string build_what(const std::string& msg);
};

//! Exception for missing files.
class FileNotFoundException: public Exception
{
    public:
        FileNotFoundException(std::string filename);
};

//! Exception for out-of-bounds interpolation.
class InterpOutOfRangeException: public Exception
{
    public:
        InterpOutOfRangeException(double x);
};

//! Exception for wrong # of command-line arguments.
class WrongCLIUsageException: public Exception
{
    public:
        WrongCLIUsageException();
};

#endif
