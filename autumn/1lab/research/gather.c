#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MSG_SIZE 1024 * 1024

double make_gather(int rank, int commsize, char* sbuf, char* rbuf, int size)
{
    MPI_Status status;

    double t = MPI_Wtime();

    if (rank == 0)
        strncpy(rbuf, sbuf, size);

    if (rank > 0)
        MPI_Send(sbuf, size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    if (rank == 0)
        for (int i = 1; i < commsize; i++)
            MPI_Recv(rbuf + size * i, size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);

    t = MPI_Wtime() - t;
    double tmax = 0;
    MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    return tmax;
}

int main(int argc, char** argv)
{
    int size = MSG_SIZE;
    int rank, commsize;
    char *sbuf = malloc(size), *rbuf = NULL;

    if (!sbuf)
        return 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    char filename[1024];
    sprintf(filename, "gather_result_%d_proc.txt", commsize);
    FILE* file;
    if (rank == 0)
        file = fopen(filename, "a");

    if (rank == 0) {
        rbuf = malloc(size * commsize);
        if (!rbuf)
            return 1;
    }
    double t;

    t = make_gather(rank, commsize, sbuf, rbuf, 1024);
    if (rank == 0)
        fprintf(file, "%-10d%f\n", 1024, t);
    t = make_gather(rank, commsize, sbuf, rbuf, 1024 * 1024);
    if (rank == 0)
        fprintf(file, "%-10d%f\n\n", 1024 * 1024, t);

    if (rank == 0)
        fclose(file);

    MPI_Finalize();

    free(sbuf);
    if (rank == 0)
        free(rbuf);
}