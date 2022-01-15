// This script reads HDF5 files for use in the coherent beamformer thread
// To compile it:
// gcc read_hdf5_file.c -o read_hdf5_file.exe -lm -lhdf5
// To run it:
// ./read_hdf5_file.exe

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include <stdint.h>
#include <endian.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include "hdf5.h"

#define FILE "/datag/users/mruzinda/hdf5/blk4.uvh5"

int main()
{
/*
  hid_t file_id, dataset_id; // identifiers //
  herr_t status;

  double dset_data[64];
  // Open an existing file. //
  file_id = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open an existing dataset. //
  dataset_id = H5Dopen2(file_id, "/Header/freq_array", H5P_DEFAULT);

  // Read the dataset. //
  status = H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);

  printf("dset_data[0] = %lf \n",dset_data[0]);
  printf("dset_data[1] = %lf \n",dset_data[1]);

  // Close the dataset. //
  status = H5Dclose(dataset_id);

  // Close the file. //
  status = H5Fclose(file_id);
*/

/*
  hid_t file_id, dataset_id; // identifiers //
  herr_t status;

  double dset_data[64][3];
  // Open an existing file. //
  file_id = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open an existing dataset. //
  dataset_id = H5Dopen2(file_id, "/Header/antenna_positions", H5P_DEFAULT);

  // Read the dataset. //
  status = H5Dread(dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);

  printf("dset_data[0][0] = %lf \n",dset_data[0][0]);
  printf("dset_data[0][1] = %lf \n",dset_data[0][1]);

  // Close the dataset. //
  status = H5Dclose(dataset_id);

  // Close the file. //
  status = H5Fclose(file_id);
*/

  hid_t file_id, dataset_id; // identifiers //
  herr_t status;

  typedef struct complex_t{
    int32_t re;
    int32_t im;
  }complex_t;

  // complex_t dset_data[17110][64][4];
  complex_t *dset_data;
  dset_data = malloc(17110*64*4*sizeof(complex_t));

  hid_t reim_tid;
  reim_tid = H5Tcreate(H5T_COMPOUND, sizeof(complex_t));
  H5Tinsert(reim_tid, "r", HOFFSET(complex_t, re), H5T_STD_I32LE);
  H5Tinsert(reim_tid, "i", HOFFSET(complex_t, im), H5T_STD_I32LE);

  // Open an existing file. //
  file_id = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

  // Open an existing dataset. //
  dataset_id = H5Dopen(file_id, "/Data/visdata", H5P_DEFAULT);

  // Read the dataset. //
  status = H5Dread(dataset_id, reim_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
  //status = H5Dread(dataset_id, H5T_COMPOUND, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
  printf("dset_data[10].re = %d \n",dset_data[10].re);

  // Close the dataset. //
  status = H5Dclose(dataset_id);

  // Close the file. //
  status = H5Fclose(file_id);

  return 0;
}
