#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define MSG_SIZE 1024 * 1024

int main(int argc, char** argv)
{
    int size = MSG_SIZE;
    int rank, commsize;
    char *sbuf = malloc(size), *rbuf;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    sbuf[0] = '0' + rank;

    MPI_Request sr, rr;

    if (rank == 0)
        rbuf = malloc(size * commsize);
    MPI_Isend(sbuf, size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &sr);
    if (rank == 0)
        for (int i = 0; i < commsize; i++)
            MPI_Irecv(rbuf + size * i, size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &rr);

    MPI_Finalize();
    if (rank == 0)
        for (int i = 0; i < commsize; i++)
            printf("%d rank got message: %c\n", rank, (rbuf + size * i)[0]);
}