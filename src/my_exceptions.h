#include <string>
#include <exception>
#include <vector>

using namespace std;

//! Exception for missing files.
class FileNotFoundException: public exception
{
    public:
        FileNotFoundException(string filename);
        ~FileNotFoundException() throw();
        virtual const char* what() const throw();
    private:
        string errmsg; //!< Error message.
};

//! Exception for out-of-bounds interpolation.
class InterpOutOfRangeException: public exception
{
    public:
        InterpOutOfRangeException(vector< pair<double, double> > table, double x);
        ~InterpOutOfRangeException() throw();
        virtual const char* what() const throw();
    private:
        string errmsg; //!< Error message.
};
