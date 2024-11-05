#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

enum { m = 45000, n = 45000 };

void* xmalloc(size_t n)
{
    void* ptr = malloc(n);

    if (ptr == NULL) {
        fprintf(stderr, "Memory allocation error\n");
        exit(1);
    }

    return ptr;
}

void get_chunk(int a, int b, int commsize, int rank, int* lb, int* ub)
{
    int n = b - a + 1;
    int q = n / commsize;
    if (n % commsize)
        q++;
    int r = commsize * q - n;
    int chunk = q;
    if (rank >= commsize - r)
        chunk = q - 1;
    *lb = a;
    if (rank > 0) {
        if (rank <= commsize - r)
            *lb += q * rank;
        else
            *lb += q * (commsize - r) + (q - 1) * (rank - (commsize - r));
    }
    *ub = *lb + chunk - 1;
}

void dgemv(float* a, float* b, float* c, int m, int n)
{
    int commsize, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int lb, ub;
    get_chunk(0, m - 1, commsize, rank, &lb, &ub);
    int nrows = ub - lb + 1;
    for (int i = 0; i < nrows; i++) {
        c[lb + i] = 0.0;
        for (int j = 0; j < n; j++)
            c[lb + i] += a[i * n + j] * b[j];
    }

    int* displs = malloc(sizeof(*displs) * commsize);
    int* rcounts = malloc(sizeof(*rcounts) * commsize);
    for (int i = 0; i < commsize; i++) {
        int l, u;
        get_chunk(0, m - 1, commsize, i, &l, &u);
        rcounts[i] = u - l + 1;
        displs[i] = (i > 0) ? displs[i - 1] + rcounts[i - 1] : 0;
    }
    MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_FLOAT, c, rcounts, displs, MPI_FLOAT, MPI_COMM_WORLD);
}

int main(int argc, char** argv)
{
    FILE* file;

    int commsize, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double t = MPI_Wtime();

    int lb, ub;
    get_chunk(0, m - 1, commsize, rank, &lb, &ub);
    int nrows = ub - lb + 1;
    float* a = xmalloc(sizeof(*a) * nrows * n);
    float* b = xmalloc(sizeof(*b) * n);
    float* c = xmalloc(sizeof(*c) * m);

    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < n; j++)
            a[i * n + j] = lb + i + 1;
    }
    for (int j = 0; j < n; j++)
        b[j] = j + 1;
    dgemv(a, b, c, m, n);
    t = MPI_Wtime() - t;

    double tmax = 0;

    MPI_Reduce(&t, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        file = fopen("data.txt", "a");
        fprintf(file, "%f\t%d\n", tmax, commsize);
        fclose(file);
    }
    free(a);
    free(b);
    free(c);
    MPI_Finalize();
    return 0;
}