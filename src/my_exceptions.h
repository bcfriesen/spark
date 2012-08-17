#include <string>
#include <exception>

using namespace std;

class FileNotFoundException: public exception
{
    public:
        FileNotFoundException(string filename);
        ~FileNotFoundException() throw();
        virtual const char* what() const throw();
    private:
        string errmsg;
};
