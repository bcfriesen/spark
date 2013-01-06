#include <iostream>
#include <fstream>
#include <string>
#include <yaml-cpp/yaml.h>
#include <my_exceptions.hpp>
#include <params.hpp>

using namespace std;
void read_params(string yaml_file, ParamsClass* params)
{
    /* YAML file has all the file locations, parameters, etc. */
    ifstream infile(yaml_file.c_str(), ifstream::in);

    if (!infile) throw FileNotFound(yaml_file);

    // I stole this try...catch segment straight from the YAML tutorial.
    try
    {
        YAML::Parser parser(infile);
        YAML::Node doc;
        parser.GetNextDocument(doc);
        // we read PHOENIX-style layer files
        doc["layer_file"] >> params->layer_file;
        // read # of core-intersecting rays (this is a knob)
        doc["num_core_intersect_rays"] >> params->num_core_intersect_rays;
        // get name of file for dumping output messages
        doc["output_file"] >> params->output_file;
    }
    catch(YAML::ParserException& e)
    {
        cout << e.what() << endl;
    }

    infile.close();

    return;
}
