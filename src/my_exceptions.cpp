#include <my_exceptions.h>

FileNotFoundException::FileNotFoundException(string filename)
{
    errmsg = "File not found: " + filename;
}

FileNotFoundException::~FileNotFoundException() throw() { }

const char* FileNotFoundException::what() const throw()
{
    return this->errmsg.c_str();
}
