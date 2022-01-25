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

//#define FILE "/datag/users/mruzinda/hdf5/blk4.uvh5"
#define FILE "/datag/users/mruzinda/hdf5/test.bfr5"

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

  hid_t file_id, npol_id, nbeams_id, obs_id, cal_all_id, delays_id, rates_id, time_array_id, sid1, sid2, sid3, sid4, obs_type, native_obs_type; // identifiers //
  herr_t status, cal_all_elements, delays_elements, rates_elements, time_array_elements;

  typedef struct complex_t{
    float re;
    float im;
  }complex_t;

  hid_t reim_tid;
  reim_tid = H5Tcreate(H5T_COMPOUND, sizeof(complex_t));
  H5Tinsert(reim_tid, "r", HOFFSET(complex_t, re), H5T_IEEE_F32LE);
  H5Tinsert(reim_tid, "i", HOFFSET(complex_t, im), H5T_IEEE_F32LE);

  complex_t *cal_all_data;
  double *delays_data;
  double *rates_data;
  double *time_array_data;
  uint64_t nbeams;
  uint64_t npol;

  int Nant = 61;    // Number of antennas
  int Nbeams = 61;  // Number of beams
  int Ntimes = 300; // Number of time stamps
  int Npol = 2;     // Number of polarizations
  int a = 34; // Antenna index
  int b = 1;  // Beam index
  int t = 1;  // Time stamp index
  int p = 1;  // Polarization index
  int c = 223;// Coarse channel index

  // Open an existing file. //
  file_id = H5Fopen(FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

  // -------------Read obsid first----------------- //
  // Open an existing dataset. //
  obs_id = H5Dopen(file_id, "/obsinfo/obsid", H5P_DEFAULT);
  // Get obsid data type //
  obs_type = H5Dget_type(obs_id);
  native_obs_type = H5Tget_native_type(obs_type, H5T_DIR_DEFAULT);
  int obsid_strsize = (int)H5Tget_size(native_obs_type);
  printf("obsid string size = %d\n", obsid_strsize);
  // Allocate memory to string array
  char obsid_str[obsid_strsize+1];
  obsid_str[obsid_strsize] = '\0'; // Terminate string
  // Read the dataset. //
  status = H5Dread(obs_id, native_obs_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, obsid_str);
  printf("obsid = %s \n", obsid_str);
  // Close the dataset. //
  status = H5Dclose(obs_id);
  // -----------------------------------------------//

  // Open an existing dataset. //
  cal_all_id = H5Dopen(file_id, "/calinfo/cal_all", H5P_DEFAULT);
  delays_id = H5Dopen(file_id, "/delayinfo/delays", H5P_DEFAULT);
  rates_id = H5Dopen(file_id, "/delayinfo/rates", H5P_DEFAULT);
  time_array_id = H5Dopen(file_id, "/delayinfo/time_array", H5P_DEFAULT);
  npol_id = H5Dopen(file_id, "/diminfo/npol", H5P_DEFAULT);
  nbeams_id = H5Dopen(file_id, "/diminfo/nbeams", H5P_DEFAULT);

  // Get dataspace ID //
  sid1 = H5Dget_space(cal_all_id);
  sid2 = H5Dget_space(delays_id);
  sid3 = H5Dget_space(rates_id);
  sid4 = H5Dget_space(time_array_id);

  // Gets the number of elements in the data set //
  cal_all_elements=H5Sget_simple_extent_npoints(sid1);
  delays_elements=H5Sget_simple_extent_npoints(sid2);
  rates_elements=H5Sget_simple_extent_npoints(sid3);
  time_array_elements=H5Sget_simple_extent_npoints(sid4);
  printf("Number of elements in the cal_all dataset is : %d\n", cal_all_elements);
  printf("Number of elements in the delays dataset is : %d\n", delays_elements);
  printf("Number of elements in the rates dataset is : %d\n", rates_elements);
  printf("Number of elements in the time_array dataset is : %d\n", time_array_elements);

  // Allocate memory for array
  cal_all_data = malloc((int)cal_all_elements*sizeof(complex_t));
  delays_data = malloc((int)delays_elements*sizeof(double));
  rates_data = malloc((int)rates_elements*sizeof(double));
  time_array_data = malloc((int)time_array_elements*sizeof(double));

  // Read the dataset //
  status = H5Dread(cal_all_id, reim_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, cal_all_data);
  printf("cal_all_data[%d].re = %f \n", a + Nant*p + Npol*Nant*c, cal_all_data[a + Nant*p + Npol*Nant*c].re);

  status = H5Dread(delays_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, delays_data);
  printf("delays_data[%d] = %lf \n", a + Nant*b + Nbeams*Nant*t, delays_data[a + Nant*b + Nbeams*Nant*t]);

  status = H5Dread(rates_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rates_data);
  printf("rates_data[%d] = %lf \n", a + Nant*b + Nbeams*Nant*t, rates_data[a + Nant*b + Nbeams*Nant*t]);

  status = H5Dread(time_array_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_array_data);
  printf("time_array_data[0] = %lf \n", time_array_data[0]);

  status = H5Dread(npol_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npol);
  printf("npol = %d \n", (int)npol);

  status = H5Dread(nbeams_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nbeams);
  printf("nbeams = %d \n", (int)nbeams);

  // Close the dataset. //
  status = H5Dclose(cal_all_id);
  status = H5Dclose(delays_id);
  status = H5Dclose(rates_id);
  status = H5Dclose(time_array_id);
  status = H5Dclose(npol_id);
  status = H5Dclose(nbeams_id);

  // Close the file. //
  status = H5Fclose(file_id);

  // Ordering/indexing of the 3D array //
  int aa[2][2][2];
  for(int i = 0; i < 2; i++){
    for(int j = 0; j < 2; j++){
      for(int k = 0; k < 2; k++){
        aa[i][j][k] = k + 2*j + 2*2*i;
        printf("a[%d][%d][%d] = %d\n", i, j, k, aa[i][j][k]);
      }
    }
  }

  return 0;
}
