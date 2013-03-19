#include <iostream>
#include <fstream>
#include <string>
#include <yaml-cpp/yaml.h>
#include <my_exceptions.hpp>
#include <params.hpp>

using namespace std;
void read_params(string yaml_file, ParamsClass* params)
{
    try
    {
        YAML::Node doc(yaml_file.c_str());
        // we read PHOENIX-style layer files
        params->layer_file = doc["layer_file"].as<std::string>();
        // read # of core-intersecting rays (this is a knob)
        params->num_core_intersect_rays = doc["num_core_intersect_rays"].as<unsigned int>();
        // get name of file for dumping output messages
        params->output_file = doc["output_file"].as<std::string>();
    }
    catch(YAML::ParserException& e)
    {
        cout << e.what() << endl;
    }

    return;
}
