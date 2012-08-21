#include <iostream>
#include <grid.h>
#include <misc.h>
#include <my_exceptions.h>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 2) throw WrongCLIUsageException();

    GridClass grid(argv[1]);;

    return 0;
}
