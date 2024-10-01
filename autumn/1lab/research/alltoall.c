#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MSG_SIZE 1024 * 1024

double make_alltoall(int rank, int commsize, char* sbuf, char* rbuf, int size)
{
    MPI_Request sendr[commsize - 1], recvr[commsize - 1];

    double t = MPI_Wtime();

    strncpy(rbuf + size * rank, sbuf, size);
    int j = 0;

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

    return tmax;
}

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

    if (!sbuf || !rbuf)
        return 1;

    char filename[1024];
    sprintf(filename, "alltoall_result_%d_proc.txt", commsize);
    FILE* file;
    if (rank == 0)
        file = fopen(filename, "a");

    double t;

    t = make_alltoall(rank, commsize, sbuf, rbuf, 1024);
    if (rank == 0)
        fprintf(file, "%-10d%f\n", 1024, t);
    t = make_alltoall(rank, commsize, sbuf, rbuf, 1024 * 1024);
    if (rank == 0)
        fprintf(file, "%-10d%f\n\n", 1024 * 1024, t);

    if (rank == 0)
        fclose(file);

    MPI_Finalize();

    free(sbuf);
    free(rbuf);
}