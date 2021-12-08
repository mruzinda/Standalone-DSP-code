// C program to implement read from FIFO
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#define N_DELAYS 8192
#define N_ANT 64
//#define N_BEAM 64
#define DELAY_POLYS 2
#define delay_idx(d, a, b)  (d + DELAY_POLYS*a + DELAY_POLYS*N_ANT*b) // Should be correct indexing
//#define delay_idx(d, a, b)  (b + N_BEAM*a + N_BEAM*N_ANT*d) // Should be correct indexing

int main()
{
    // File descriptor
    int fd1;

    // Array of floats to place data read from file
    float delay_pols[N_DELAYS];

    // FIFO file path
    char * myfifo = "/tmp/katpoint_delays";

    // Creating the named file(FIFO)
    // mkfifo(<pathname>,<permission>)
    mkfifo(myfifo, 0666);

    // Open file as a read only
    fd1 = open(myfifo,O_RDONLY);

    // Read file
    read(fd1, delay_pols, sizeof(delay_pols));

    // Close file
    close(fd1);

    printf("Size of delay array %lu\n", sizeof(delay_pols));

    // First beam
    printf("--------------First beam delay offset---------------\n");
    printf("idx %d in result array = %e \n", delay_idx(0, 0, 0), delay_pols[delay_idx(0, 0, 0)]);
    printf("idx %d in result array = %e \n", delay_idx(0, 1, 0), delay_pols[delay_idx(0, 1, 0)]);
    printf("idx %d in result array = %e \n", delay_idx(0, 2, 0), delay_pols[delay_idx(0, 2, 0)]);
    // Second beam delay
    printf("--------------Second beam delay offset--------------\n");
    printf("idx %d in result array = %e \n", delay_idx(0, 0, 1), delay_pols[delay_idx(0, 0, 1)]);
    printf("idx %d in result array = %e \n", delay_idx(0, 1, 1), delay_pols[delay_idx(0, 1, 1)]);
    printf("idx %d in result array = %e \n", delay_idx(0, 2, 1), delay_pols[delay_idx(0, 2, 1)]);
    // Second beam rate
    printf("---------------Second beam delay rate----------------\n");
    printf("idx %d in result array = %e \n", delay_idx(1, 0, 1), delay_pols[delay_idx(1, 0, 1)]); // 129
    printf("idx %d in result array = %e \n", delay_idx(1, 1, 1), delay_pols[delay_idx(1, 1, 1)]); // 131
    printf("idx %d in result array = %e \n", delay_idx(1, 2, 1), delay_pols[delay_idx(1, 2, 1)]); // 133

    return 0;
}
