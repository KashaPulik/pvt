#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define MSG_SIZE 1024 * 1024

double make_broadcast(int rank, int commsize, char* sbuf, char* rbuf, int size)
{
    MPI_Request sr, rr;

    double t = MPI_Wtime();

    for (int i = 0; i < commsize; i++) {
        if (rank == 0)
            MPI_Isend(sbuf, size, MPI_CHAR, i, 0, MPI_COMM_WORLD, &sr);
    }
    MPI_Irecv(rbuf, size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &rr);

    t = MPI_Wtime() - t;
    double tmax = 0;
    MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    return tmax;
}

int main(int argc, char** argv)
{
    int size = MSG_SIZE;
    int rank, commsize;
    char* sbuf = malloc(size);
    char* rbuf = malloc(size);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    char filename[1024];
    sprintf(filename, "broadcast_result_%d_proc.txt", commsize);
    FILE* file;
    if (rank == 0)
        file = fopen(filename, "a");

    double t;

    t = make_broadcast(rank, commsize, sbuf, rbuf, 1024);
    if (rank == 0)
        fprintf(file, "%-10d%f\n", 1024, t);
    t = make_broadcast(rank, commsize, sbuf, rbuf, 1024 * 1024);
    if (rank == 0)
        fprintf(file, "%-10d%f\n\n", 1024 * 1024, t);

    if (rank == 0)
        fclose(file);

    MPI_Finalize();
    free(sbuf);
    free(rbuf);
}