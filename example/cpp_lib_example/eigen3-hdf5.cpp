#include <eigen3-hdf5.hpp>
#include <Eigen/Dense>
#include <hdf5.h>
#include <H5Cpp.h>
#include <iostream>
using namespace std;
void save_matrix()
{
    Eigen::Matrix3d mat;
    mat << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    H5::H5File file("filename1.h5", H5F_ACC_TRUNC);
    char file_name[];
    sprintf(file_name, "MatrixDataSetName%d",1);
    EigenHDF5::save(file, file_name, mat);
    EigenHDF5::save(file, "MatrixDataSetName", mat);
}

void load_vector()
{
    Eigen::Matrix3d mat, mat2;
    H5::H5File file("filename1.h5", H5F_ACC_RDONLY);
    char file_name[];
    sprintf(file_name, "MatrixDataSetName%d",1);
    EigenHDF5::load(file, file_name, mat);
    EigenHDF5::load(file, "MatrixDataSetName", mat2);
    cout << mat << endl;
    cout << mat2 << endl;
}

int main(){
    save_matrix();
    load_vector();
    return 0;
}

