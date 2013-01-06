#ifndef READ_PARAMS_H
#define READ_PARAMS_H

#include <string>
#include <params.hpp>

/** Reads in all runtime parameters from supplied YAML file. */
void read_params(std::string yaml_file, ParamsClass* params);

#endif // READ_PARAMS_H
