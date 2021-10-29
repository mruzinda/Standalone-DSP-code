// This script reads the filterbank files gets a block of data and writes it to a text file for further processing.
// To compile it:
// gcc fb_file_analysis.c -o fb_file_analysis.exe -I/home/mruzinda/beamformer_workspace/include -L/home/mruzinda/beamformer_workspace/lib/ -lrawspec
// To run it:
// ./fb_file_analysis.exe /datag/users/mruzinda/oics/guppi_59143_55142_000486_GPS-BIIR-11_0001-ics.rawspec.0000.fil
// where the file in the argument is an example
// Remember to change the ics_flag in the code to 1 if you are processing incoherent sum data.
// And change the filename of the text file if necessary.

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

int main()
{
	// 1 block is 3 floats with the first being the "header", and the next 2 being the payload
	int full_size = 15;
	int blk_size = 3;
	int num_blks = full_size/blk_size;
	int hdr_size = 1;
	int payload_size = blk_size-hdr_size;
	//int only_payload = num_blks*payload_size;
	int count = 0;

	char raw_filename[128];
	char txt_filename[128];

	float * buff  = (float *)calloc(full_size, sizeof(float));
	float * buff2 = (float *)calloc(full_size, sizeof(float));
	float * buff3 = (float *)calloc(full_size, sizeof(float));

	for(int i = 0; i < full_size; i++){
		count += 1;
		buff[i] = count;
		printf("idx %d in buff val = %g\n",i, buff[i]);
	}

	strcpy(raw_filename, "raw_test.bin"); 
	strcpy(txt_filename, "txt_test.txt");

	int raw_file;
	raw_file = open(raw_filename, O_RDWR|O_CREAT);

	// Write to raw file
	printf("=========== Write to raw file ===================\n");
	write(raw_file, buff, full_size*sizeof(float));

	close(raw_file);
	free(buff);

	// Read from raw file
	printf("=========== Read from raw file ===================\n");
	raw_file = open(raw_filename, O_RDONLY);
	/*
	read(raw_file, buff2, full_size*sizeof(float));
	for(int i = 0; i < full_size; i++){
		printf("idx %d in buff val = %g\n",i, buff2[i]);
	}
	*/
	// Read 1 block at a time - Read the "header" first, then the payload
	for(int i = 0; i < num_blks; i++){
		read(raw_file, &buff2[i*blk_size], hdr_size*sizeof(float));
		read(raw_file, &buff2[(i*blk_size)+hdr_size], payload_size*sizeof(float));
	}
	for(int i = 0; i < full_size; i++){
		printf("idx %d in buff val = %g\n",i, buff2[i]);
	}
	close(raw_file);

	// Write to text file
	printf("=========== Write to text file ===================\n");
	FILE* txt_file;
	txt_file = fopen(txt_filename, "w");
	for(int i = 0; i < full_size; i++){
		printf("idx %d in buff val = %g\n",i, buff2[i]);
		fprintf(txt_file, "%g\n", buff2[i]);
	}
	fclose(txt_file);
	free(buff2);

	// Read from text file
	printf("=========== Read from text file ===================\n");
	txt_file = fopen(txt_filename, "r");
	for(int i = 0; i < full_size; i++){
		fscanf(txt_file, "%g\n", &buff3[i]);
		printf("idx %d in buff val = %g\n",i, buff3[i]);
	}
	fclose(txt_file);
	free(buff3);

	return 0;
}
