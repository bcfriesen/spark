#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <H5Cpp.h>
#include <boost/lexical_cast.hpp>
#include <my_exceptions.h>
#include <make_linelist_hdf.h>

using namespace std;

int make_linelist_hdf(const char* filename)
{
    const H5std_string HDF_FILE_NAME("lines.h5");
    const H5std_string HDF_DATASET_NAME("KuruczLines");

    /* Dimension of the data space. Ours is a vector of quantities related to a
     * single line, so RANK = 1. */
    const int HDF_RANK = 1;

    /* Number of components in the data structure. */
    const int HDF_LENGTH = 28;
    hsize_t line_dim[1];
    line_dim[0] = HDF_LENGTH;

    /* Create 10-character-long string for term labels */
    H5::StrType strtype10(H5::PredType::C_S1, 10);
    /* Create 4-character-long strong for line data reference */
    H5::StrType strtype4(H5::PredType::C_S1, 4);

    /* Create the HDF file. */
    H5::H5File hdf_file(HDF_FILE_NAME, H5F_ACC_TRUNC);

    /* Create the data space. */
    H5::DataSpace hdf_dataspace(HDF_RANK, line_dim);

    /* Create the HDF data type. */
    H5::CompType hdf_line(sizeof(LineStruct));

    /* Fill in data type. */
    hdf_line.insertMember("Wavelength (nm)", HOFFSET(LineStruct, wl), H5::PredType::NATIVE_DOUBLE);
    hdf_line.insertMember("Log10(gf) (cgs)", HOFFSET(LineStruct, loggf), H5::PredType::NATIVE_DOUBLE);
    hdf_line.insertMember("Kurucz element code", HOFFSET(LineStruct, elcode), H5::PredType::NATIVE_DOUBLE);
    hdf_line.insertMember("Upper level energy (cm^-1)", HOFFSET(LineStruct, energy1), H5::PredType::NATIVE_DOUBLE);
    hdf_line.insertMember("Lower level energy (cm^-1)", HOFFSET(LineStruct, energy2), H5::PredType::NATIVE_DOUBLE);
    hdf_line.insertMember("Upper level angular momentum", HOFFSET(LineStruct, J1), H5::PredType::NATIVE_DOUBLE);
    hdf_line.insertMember("Lower level angular momentum", HOFFSET(LineStruct, J2), H5::PredType::NATIVE_DOUBLE);

    /* Create the dataset. */
    H5::DataSet hdf_dataset = hdf_file.createDataSet(HDF_DATASET_NAME, hdf_line, hdf_dataspace);


    vector<LineStruct> lines;
    LineStruct line;

    char blank[10];
    string oneline;

    /* Read Kurucz line data. */
    ifstream infile;
    infile.open(filename);
    if (!infile) throw FileNotFound(filename);

    while(!infile.eof())
    {
        getline(infile, oneline);
        istringstream is(oneline);

        is >> line.wl >> line.loggf >> line.elcode >> line.energy1 >> line.J1 >> blank >> line.energy2 >> line.J2;

        cout << "wl = " << line.wl << endl;
        cout << "loggf = " << line.loggf << endl;
        cout << "elcode = " << line.elcode << endl;
        cout << "energy1 = " << line.energy1 << endl;
        cout << "energy2 = " << line.energy2 << endl;
        cout << "J1 = " << line.J1 << endl;
        cout << "J2 = " << line.J2 << endl;
        cout << endl;

        lines.push_back(line);

        cout << "writing to h5 file..." << endl;
        hdf_dataset.write(&line, hdf_line);
    }

    infile.close();

    return 0;
}
