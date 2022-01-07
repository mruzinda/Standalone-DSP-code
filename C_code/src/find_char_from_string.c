// This script tests reading and writing blocks of data to and from both RAW and text files
// Also using this script to test python module importing.
// To compile it:
// gcc find_char_from_string.c -o find_char_from_string.exe -lm
// To run it:
// ./find_char_from_string.exe

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

#define ARRAY_SIZE 100

int main()
{
  char *fname;
  char path[ARRAY_SIZE] = {0};
  fname = "/datag/users/mruzinda/i/guppi_59143_55142_000486_GPS-BIIR-11_0001.0000.raw";
  printf("File name: %s\n",fname);

  // Find path to file in file name
  char character = '/';

  // strrchr() finds the last occurence of the specified character
  char *ptr; 
  long int char_pos;
  ptr = strrchr(fname, character);
  char_pos = ptr-fname;
  printf("The last position of %c is %ld \n", character, char_pos);

  // Copy path portion of file name to new variable
  memcpy(path, fname, char_pos+1);

  printf("Path to files is: %s\n", path);

  // Go to the directory and look for a file with "0000.raw" and wait until it shows up
  char first_raw_ext[ARRAY_SIZE] = {0};
  char first_raw[ARRAY_SIZE] = {0};
  char command[ARRAY_SIZE] = {0};
  strcpy(first_raw_ext, "*0000.raw");
  strcpy(first_raw, path);
  strcat(first_raw, first_raw_ext);
  strcpy(command, "ls ");
  strcat(command, first_raw);

  printf("Ubuntu command: %s\n", command);

  FILE *fp;
  int status;
  char filename[ARRAY_SIZE] = {0};


  fp = popen(command, "r");
  if (fp == NULL){
    printf("Failed to execute command!\n");
  }

  if (fgets(filename, ARRAY_SIZE, fp) != NULL){
    if(strncmp(fname,filename,strlen(fname)) == 0){
      printf("Got filename: %s\n", filename);
    }
  }
  else{
    printf("There's nothing in this directory! \n");
  }

  status = pclose(fp);
  if (status == -1) {
    printf("Failed to close pipe!\n");
  }

  return 0;
}
