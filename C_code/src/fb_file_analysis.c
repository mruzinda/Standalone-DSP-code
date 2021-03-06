// This script reads the filterbank files gets a block of data and writes it to a text file for further processing.
// To compile it:
// gcc fb_file_analysis.c -o fb_file_analysis.exe -I/home/mruzinda/hpguppi_proc/include -L/home/mruzinda/hpguppu_proc/lib/ -lrawspec
// To run it:
// ./fb_file_analysis.exe /datag/users/mruzinda/oics/guppi_59143_55142_000486_GPS-BIIR-11_0001-ics.rawspec.0000.fil
// where the file in the argument is an example
// Remember to change the sinlge_beam_flag in the code to 1 if you are processing incoherent sum data.
// And change the filename of the text file if necessary.

#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <endian.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include "rawspec_fbutils.h"
#include "coherent_beamformer_char_in.h"

int main(int argc, char * argv[])
{
  	fb_hdr_t hdr;

	// 0=fake data; 1=Arecibo; 2=Ooty... others to be added
  	hdr.machine_id = 20; // wtf?
  	// 0=FAKE; 1=PSPM; 2=WAPP; 3=OOTY... others to be added
  	hdr.telescope_id = 6; // GBT
  	// 1=filterbank; 2=time series... others to be added
  	hdr.data_type = 1;
  	// 1 if barycentric or 0 otherwise (only output if non-zero)
  	hdr.barycentric = 1;
  	// 1 if pulsarcentric or 0 otherwise (only output if non-zero)
  	hdr.pulsarcentric = 1;
  	// right ascension (J2000) of source (hours)
  	// will be converted to/from hhmmss.s
  	hdr.src_raj = 20.0 + 39/60.0 + 7.4/3600.0;
  	// declination (J2000) of source (degrees)
  	// will be converted to/from ddmmss.s
  	hdr.src_dej = 42.0 + 24/60.0 + 24.5/3600.0;
  	// telescope azimuth at start of scan (degrees)
  	hdr.az_start = 12.3456;
  	// telescope zenith angle at start of scan (degrees)
  	hdr.za_start = 65.4321;
  	// centre frequency (MHz) of first filterbank channel
  	hdr.fch1 = 4626.464842353016138;
  	// filterbank channel bandwidth (MHz)
  	hdr.foff = -0.000002793967724;
  	// number of filterbank channels
  	hdr.nchans = 1;
  	// total number of beams
  	hdr.nbeams = 1;
  	// total number of beams
  	hdr.ibeam = 1;
  	// number of bits per time sample
  	hdr.nbits = 32;
  	// time stamp (MJD) of first sample
  	hdr.tstart = 57856.810798611114;
  	// time interval between samples (s)
  	hdr.tsamp = 1.825361100800;
  	// number of seperate IF channels
  	hdr.nifs = 1;
  	// the name of the source being observed by the telescope
  	// Max string size is supposed to be 80, but bug in sigproc if over 79
  	strcpy(hdr.source_name, "1234567890123456789012345678901234567890123456789012345678901234567890123456789");
  	// the name of the original data file
  	// Max string size is supposed to be 80, but bug in sigproc if over 79
  	strcpy(hdr.rawdatafile, "1234567890123456789012345678901234567890123456789012345678901234567890123456789");

  	if(argc > 1) {
    		int fd  = open(argv[1],  O_RDONLY);
		ssize_t sz = 1;
		int sinlge_beam_flag = 1; // Set to 1 if the incoherent sum processed the GUPPI RAW file
		int blk_size = 0;
		int count_blks = 0; // Used to count the number of blocks in filterbank file
		if(sinlge_beam_flag == 1){
			blk_size = N_BF_POW/N_BEAM;
		}else{
			blk_size = N_BF_POW;
		}

		// Write data to text file for analysis
		char output_filename[128];

		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test.txt"); // cbf
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test1.txt"); // ics upchannelized
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test2.txt"); // ics no upchannelization
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test3.txt"); // cbf
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test4.txt"); // sim-cbf
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test5.txt"); // cbf-after simulated data was correct
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test6.txt"); // ics no upchannelization -after simulated data was correct
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test7.txt"); // cbf - RAW FILE input
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test8.txt"); // ics upchannelization - RAW FILE input
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test9.txt"); // cbf - RAW file input after working moving SOI simulated
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test10.txt"); // cbf - RAW file input after changing read_fully() back
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test11.txt"); // cbf - RAW file input after change from 32 to 64 coarse channs
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test12.txt"); // sim-cbf - After change from 32 to 64 coarse channs
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test13.txt"); // cbf - RAW file input after time averaging added to python post-proc
		//strcpy(output_filename, "/datag/users/mruzinda/out_txt/output_d_test14.txt"); // cbf - RAW file input after time averaging added to python post-proc
		strcpy(output_filename, "/mydatag/20220120/0024/output_d_test.txt"); // cbf - RAW FILE input - Only one beam per file (so same structure as rawspec)

		FILE* output_file;

		output_file = fopen(output_filename, "a"); // Append rather than overwrite file

		//int blk_plus_hdr_elem = hdr_size+blk_size;
		//int blk_plus_hdr_bytes = hdr_size+(blk_size*sizeof(float));
		//float * buff = (float *)calloc(blk_plus_hdr, sizeof(float));
		//sz = read(fd, buff, blk_plus_hdr*sizeof(float));
		//float * data_payload = (float *)(buff + hdr_size);
		float * buff = (float *)calloc(blk_size, sizeof(float));
		float * data_payload;
		ssize_t hdr_size;
		while(sz>0){
			if(sinlge_beam_flag == 1 && count_blks == 0){
				hdr_size = fb_fd_read_header(fd, &hdr, NULL);
				printf("header size %lu bytes\n", hdr_size);
    				printf("fch1 %.17g\n", hdr.fch1);
				printf("foff %.17g\n", hdr.foff);
				printf("nbeams %d\n", hdr.nbeams);
                                printf("src_raj %g\n", hdr.src_raj);
				printf("src_dej %g\n", hdr.src_dej);
			}
			if(sinlge_beam_flag == 0){ // Currently, coherent beamformer 
				hdr_size = fb_fd_read_header(fd, &hdr, NULL);
				printf("header size %lu bytes\n", hdr_size);
    				printf("fch1 %.17g\n", hdr.fch1);
				printf("foff %.17g\n", hdr.foff);
				printf("nbeams %d\n", hdr.nbeams);
                                printf("src_raj %g\n", hdr.src_raj);
				printf("src_dej %g\n", hdr.src_dej);
			}
    			
			sz = read(fd, buff, blk_size*sizeof(float));
			printf("Size of payload read from filterbank file in bytes: %zd \n", sz);
			if(sz < (blk_size*sizeof(float))) {
				printf("All bytes seem to have been read!\n");
                		break;
        		}

			data_payload = (float *)(buff);

			for (int ii = 0; ii < blk_size; ii++) {
				//fprintf(output_file, "%c\n", output_data[ii]);
				fprintf(output_file, "%g\n", data_payload[ii]);
			}
			memset(buff, 0, blk_size*sizeof(float)); // Set to zero after block is written to file (Probably not necessary)
			count_blks += 1;
			printf("Block index read and written: %d \n", count_blks);
		}
		printf("Number of blocks in the filter bank file = %d \n", count_blks);

		fclose(output_file);
		//close(argv[1]);
		printf("Closed output file.\n");

  	} else {
		printf("Must add filterbank file as argument.\n");
/*
    int i;
    float f = 0;
    char fname[80];

    for(i=0; i<100; i++) {
      sprintf(fname, "fbutils_fd.%02d.fil", i);
      int fdfd  = open(fname,  O_WRONLY | O_CREAT, 0664);
      sprintf(fname, "fbutils_buf.%02d.fil", i);
      int fdbuf = open(fname, O_WRONLY | O_CREAT, 0664);


      ssize_t nbytes = fb_fd_write_padded_header(fdfd, &hdr, 1024+i);
      write(fdfd, (void *)&f, sizeof(float));
      printf("%02d: write %lu+4 fd bytes, ", i, nbytes);

      char buf[1024];
      char * end = fb_buf_write_padded_header(buf, &hdr, 1024+i);
      nbytes = end-buf;
      write(fdbuf, buf, nbytes);
      write(fdbuf, (void *)&f, sizeof(float));
      printf("write %lu+4 buf bytes\n", nbytes);

      close(fdfd);
      close(fdbuf);
    }
*/
  	}

  return 0;
}
