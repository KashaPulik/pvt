#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define MSG_SIZE 1024 * 1024

double make_ring(int rank, int commsize, char* msg_send, char* msg_recv, int size)
{
    MPI_Status status;
    MPI_Request send_request, recv_request;

    double t = MPI_Wtime();

    for (int step = 0; step < commsize - 1; step++) {
        int next = (rank + 1) % commsize;
        int prev = (rank - 1 + commsize) % commsize;

        MPI_Isend(msg_send, size, MPI_CHAR, next, 0, MPI_COMM_WORLD, &send_request);
        MPI_Irecv(msg_recv, size, MPI_CHAR, prev, 0, MPI_COMM_WORLD, &recv_request);

        MPI_Wait(&send_request, &status);
        MPI_Wait(&recv_request, &status);
    }

    t = MPI_Wtime() - t;
    double tmax = 0;
    MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    return tmax;
}

int main(int argc, char** argv)
{
    int rank, commsize;
    char* msg_send = malloc(MSG_SIZE);
    char* msg_recv = malloc(MSG_SIZE);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    char filename[1024];
    sprintf(filename, "ring_result_%d_proc.txt", commsize);
    FILE* file;
    if (rank == 0)
        file = fopen(filename, "a");

    double t;

    t = make_ring(rank, commsize, msg_send, msg_recv, 1);
    if (rank == 0)
        fprintf(file, "%-10d%f\n", 1, t);
    t = make_ring(rank, commsize, msg_send, msg_recv, 1024);
    if (rank == 0)
        fprintf(file, "%-10d%f\n", 1024, t);
    t = make_ring(rank, commsize, msg_send, msg_recv, 1024 * 1024);
    if (rank == 0)
        fprintf(file, "%-10d%f\n\n", 1024 * 1024, t);

    if(rank == 0)
        fclose(file);
    
    MPI_Finalize();

    free(msg_send);
    free(msg_recv);
    return 0;
}
