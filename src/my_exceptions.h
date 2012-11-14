#ifndef MY_EXCEPTIONS_H
#define MY_EXCEPTIONS_H

#include <string>
#include <stdexcept>

/** \brief Error message namespace.
 *
 * Shamelessly stolen from yaml-cpp.
 */
namespace ErrorMsg
{
    const std::string FILE_NOT_FOUND = "File not found: ";
    const std::string INTERP_OUT_OF_RANGE = "Interpolation out of range for value: ";
    const std::string WRONG_CLI_USAGE = "Usage: <binary> <yaml_file>";
    const std::string NONMONOTONIC = "Non-monotonic velocity field! Error on this velocity: ";
}

/** \brief Generic exception class.
 *
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

/** Missing files. */
class FileNotFound: public Exception
{
    public:
        FileNotFound(std::string filename);
};

/** Out-of-bounds interpolation. */
class InterpOutOfRange: public Exception
{
    public:
        InterpOutOfRange(double x);
};

/** Wrong # of command-line arguments. */
class WrongCLIUsage: public Exception
{
    public:
        WrongCLIUsage();
};

/** Non-monotonic velocity field. */
class NonmonotonicVelocityField: public Exception
{
    public:
        NonmonotonicVelocityField(double vel);
};

#endif
