#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <H5Cpp.h>
#include <boost/lexical_cast.hpp>
#include <my_exceptions.hpp>
#include <make_linelist_hdf.hpp>

using namespace std;

void make_linelist_hdf(const char* filename)
{
    const H5std_string HDF_FILE_NAME("lines.h5");
    const H5std_string HDF_DATASET_NAME("KuruczLines");

    /* Dimension of the data space. Ours is a vector of quantities describing a
     * single atomic line, so RANK = 2. */
    const int HDF_RANK = 2;

    /* Number of components in the data structure. */
    const int HDF_LENGTH = 7;
    /* The components of the data structure. */
    double wl, loggf, elcode, energy1, energy2, J1, J2;

    /* The data space constructor expects the number of dimensions to be
     * contained in an array, even if there's only one dimension. */
    hsize_t line_dim[2], max_line_dim[2];
    line_dim[0] = 1;
    line_dim[1] = HDF_LENGTH;
    max_line_dim[0] = H5S_UNLIMITED;
    max_line_dim[1] = HDF_LENGTH;

    /* Create the HDF file. */
    H5::H5File hdf_file(HDF_FILE_NAME, H5F_ACC_TRUNC);

    /* In order to append the dataset, it must be "chunked." */
    H5::DSetCreatPropList hdf_cparms;
    hsize_t chunk_dims[2] = {1, HDF_LENGTH};
    hdf_cparms.setChunk(HDF_RANK, chunk_dims);

    /* Create the data space. A missing third argument means the maximum extent
     * of the data space is unlimited.*/
    H5::DataSpace hdf_dataspace(HDF_RANK, line_dim, max_line_dim);

    /* Create the dataset. */
    H5::DataSet hdf_dataset = hdf_file.createDataSet(HDF_DATASET_NAME,
                                                     H5::PredType::NATIVE_DOUBLE,
                                                     hdf_dataspace,
                                                     hdf_cparms);

    char blank[10];
    string oneline;
    double one_line_of_data[HDF_LENGTH];

    /* Read Kurucz line data. */
    ifstream infile;
    infile.open(filename);
    if (!infile) throw FileNotFound(filename);

    /* Set up offset sizes for each new write to the dataset. */
    hsize_t offset[2];
    offset[0] = 0;
    offset[1] = 0;

    /* We must extend the dataset each time we want to append new data to it. */
    hsize_t extend[2];
    extend[0] = 1;
    extend[1] = HDF_LENGTH;

    H5::DataSpace hdf_dataspace_new_line;

    cout << "Creating HDF linelist..." << endl;

    while(getline(infile, oneline))
    {
        /* Read line data from ASCII file. */
        istringstream is(oneline);

        is >> setw(11) >> wl
           >> setw( 7) >> loggf
           >> setw( 6) >> elcode
           >> setw(12) >> energy1
           >> setw( 5) >> J1
           >> setw( 1) >> blank
           >> setw(12) >> energy2
           >> setw( 5) >> J2;

        one_line_of_data[0] = wl;
        one_line_of_data[1] = loggf;
        one_line_of_data[2] = elcode;
        one_line_of_data[3] = energy1;
        one_line_of_data[4] = energy2;
        one_line_of_data[5] = J1;
        one_line_of_data[6] = J2;

        /* Extend dataset for each write. */
        hdf_dataset.extend(extend);

        /* Set up hyperslab to write one new line. */
        hdf_dataspace_new_line = hdf_dataset.getSpace();
        hdf_dataspace_new_line.selectHyperslab(H5S_SELECT_SET, line_dim, offset);

        /* Write hyperslab data. */
        hdf_dataset.write(one_line_of_data, H5::PredType::NATIVE_DOUBLE, hdf_dataspace, hdf_dataspace_new_line);

        /* Shift to the next hyperslab for a new write. */
        offset[0] += 1;
        extend[0] += 1;
    }
    cout << "HDF linelist construction complete!" << endl << endl;

    infile.close();
}
