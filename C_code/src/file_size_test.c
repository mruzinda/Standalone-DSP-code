// C program to get binary/raw file size
// In order to compile the script use the command:
// gcc file_size_test.c -o file_size_test.exe
// And run the script, enter:
// ./file_size_test.exe
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <unistd.h>

long int get_file_size(int fdin){
    off_t cur_pos = lseek(fdin, (size_t)0, SEEK_CUR);
    off_t file_size = lseek(fdin, (size_t)0, SEEK_END);
    lseek(fdin, cur_pos, SEEK_SET);
    return file_size;
}

int main(){
    int open_flags = O_RDONLY;
    char fname[256];
    strcpy(fname, "/datag/users/mruzinda/i/input_h_cufft.bin");
    int fdin = open(fname, open_flags, 0644);
    if (fdin==-1) {
        printf("Error opening file.");
    }

    long int file_size = get_file_size(fdin);
    printf("File size = %ld bytes\n", file_size);

    long int cur_pos = lseek(fdin, (size_t)0, SEEK_CUR);
    printf("Current position = %ld bytes\n", cur_pos);

    close(fdin);
    return 0;
}
