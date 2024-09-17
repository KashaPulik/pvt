#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv)
{
    int rank, commsize;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double *buf = malloc(sizeof(*buf) * 100);

    if (rank == 0) {
        MPI_Send(buf, 100, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    } else if (rank == 1) {
        MPI_Recv(buf, 100, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    free(buf);
    MPI_Finalize();
    return 0;
}