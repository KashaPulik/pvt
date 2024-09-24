#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MSG_SIZE 1024 * 1024

int main(int argc, char** argv)
{
    int size = MSG_SIZE;
    int rank, commsize;
    char *sbuf, *rbuf;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    sbuf = malloc(size);
    rbuf = malloc(size * commsize);

    sbuf[0] = '0' + rank;
    strncpy(rbuf + size * rank, sbuf, size);

    MPI_Request sendr[commsize - 1], recvr[commsize - 1];
    int j = 0;

    double t = MPI_Wtime();

    for (int i = 0; i < commsize; i++) {
        if (i == rank)
            continue;
        MPI_Isend(sbuf, size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &sendr[j]);
        MPI_Irecv(rbuf + size * i, size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &recvr[j]);
        j++;
    }
    MPI_Waitall(commsize - 1, sendr, NULL);
    MPI_Waitall(commsize - 1, recvr, NULL);

    t = MPI_Wtime() - t;
    double tmax = 0;
    MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    for (int i = 0; i < commsize; i++)
        printf("%c ", (rbuf + size * i)[0]);
    printf("\n");
    MPI_Finalize();

    if (rank == 0)
        printf("Elapsed time: %f seconds\n", tmax);
}