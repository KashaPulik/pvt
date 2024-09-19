#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define NELEMS(x) (sizeof(x) / sizeof((x)[0]))
#define MSG_SIZE 1024 * 1024

int main(int argc, char** argv)
{
    int rank, commsize, len, tag = 1;
    char host[MPI_MAX_PROCESSOR_NAME];
    char msg_send[MSG_SIZE];
    char msg_recv[MSG_SIZE];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Get_processor_name(host, &len);
    MPI_Status status;
    MPI_Request send_request, recv_request;

    snprintf(msg_send, NELEMS(msg_send), "Hello from process %d of %d on %s", rank, commsize, host);

    for (int step = 0; step < commsize - 1; step++) {
        int next = (rank + 1) % commsize;
        int prev = (rank - 1 + commsize) % commsize;

        MPI_Isend(msg_send, NELEMS(msg_send), MPI_CHAR, next, tag, MPI_COMM_WORLD, &send_request);
        MPI_Irecv(msg_recv, NELEMS(msg_recv), MPI_CHAR, prev, tag, MPI_COMM_WORLD, &recv_request);

        MPI_Wait(&send_request, &status);
        MPI_Wait(&recv_request, &status);

        snprintf(msg_send, NELEMS(msg_send), "%s", msg_recv);
        printf("Message received by process %d from process %d: '%s'\n", rank, prev, msg_recv);
    }

    MPI_Finalize();
    return 0;
}
