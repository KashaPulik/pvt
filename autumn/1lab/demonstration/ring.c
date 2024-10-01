#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define NELEMS(x) (sizeof(x) / sizeof((x)[0]))
#define MSG_SIZE 1024

int main(int argc, char** argv)
{
    int rank, commsize, len;
    char host[MPI_MAX_PROCESSOR_NAME];
    char* msg_send = malloc(MSG_SIZE);
    char* msg_recv = malloc(MSG_SIZE);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Get_processor_name(host, &len);
    MPI_Request send_request, recv_request;

    snprintf(msg_send, MSG_SIZE, "Hello from process %d of %d on %s", rank, commsize, host);

    for (int step = 0; step < commsize - 1; step++) {
        int next = (rank + 1) % commsize;
        int prev = (rank - 1 + commsize) % commsize;

        MPI_Isend(msg_send, MSG_SIZE, MPI_CHAR, next, 0, MPI_COMM_WORLD, &send_request);
        MPI_Irecv(msg_recv, MSG_SIZE, MPI_CHAR, prev, 0, MPI_COMM_WORLD, &recv_request);

        MPI_Waitall(1, &send_request, NULL);
        MPI_Waitall(1, &recv_request, NULL);

        snprintf(msg_send, MSG_SIZE, "%s", msg_recv);
        printf("Message received by process %d from process %d: '%s'\n", rank, prev, msg_recv);
    }

    MPI_Finalize();
    free(msg_send);
    free(msg_recv);
    return 0;
}
