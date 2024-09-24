#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MSG_SIZE 1024

int main(int argc, char** argv)
{
    int size = MSG_SIZE;
    int rank, commsize;
    char *sbuf = malloc(size), *rbuf;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    sbuf[0] = '0' + rank;

    MPI_Status status;

    double t = MPI_Wtime();
    if (rank == 0) {
        rbuf = malloc(size * commsize);
        strncpy(rbuf, sbuf, size);
    }
    if (rank > 0)
        MPI_Send(sbuf, size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    if (rank == 0)
        for (int i = 1; i < commsize; i++)
            MPI_Recv(rbuf + size * i, size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);

    t = MPI_Wtime() - t;
    double tmax = 0;
    MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    if (rank == 0)
        for (int i = 0; i < commsize; i++)
            printf("%d rank got message: %c\n", rank, (rbuf + size * i)[0]);
    free(sbuf);
    if (rank == 0)
        free(rbuf);
    if (rank == 0)
        printf("Elapsed time: %f seconds\n", tmax);
}