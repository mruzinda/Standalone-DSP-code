// C program to implement read from FIFO
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#define N_DELAYS 8192

int main()
{
    //int fd1;
    FILE *fd1;

    // FIFO file path
    char * myfifo = "/datag/users/mruzinda/katpoint_delays";

    // Creating the named file(FIFO)
    // mkfifo(<pathname>,<permission>)
    mkfifo(myfifo, 0666);

    char str1[N_DELAYS];
    
    fd1 = fopen(myfifo,"r");

    int c = getc(fd1);
    while(c != EOF)
    {
        // First open in read only and read
        //fd1 = open(myfifo,O_RDONLY);
        //read(fd1, str1, N_DELAYS);
	
	fread(str1, N_DELAYS, 1, fd1);

        // Print the read string and close
        printf("User1: %s\n", &str1[0]);
        putchar(c);
        c = getc(fd1);
        
    }
    fclose(fd1);
    return 0;
}
