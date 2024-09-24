#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define MSG_SIZE 1024 * 1024

void make_broadcast(int size, int argc, char** argv)
{
    int rank, commsize;
    char* sbuf = malloc(size);
    char* rbuf = malloc(size);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Request sr, rr;

    sbuf[0] = '0' + rank;
    for (int i = 0; i < commsize; i++) {
        if (rank == 0)
            MPI_Isend(sbuf, size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &sr);
    }
    MPI_Irecv(rbuf, size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &rr);

    MPI_Finalize();
    printf("%d rank got message: %c\n", rank, rbuf[0]);
    free(sbuf);
    free(rbuf);
}

int main(int argc, char** argv)
{
    make_broadcast(MSG_SIZE, argc, argv);
}