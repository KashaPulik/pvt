#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

const double PI = 3.14159265358979323846;
const int n = 10000000;

double getrand()
{
    return (double)rand() / RAND_MAX;
}

double func(double x, double y)
{
    return exp(x - y);
}

int main(int argc, char** argv)
{
    FILE* file;
    MPI_Init(&argc, &argv);
    int rank, commsize;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    int in = 0;
    double s = 0;
    int gin = 0;
    double gsum = 0.0;

    double time = MPI_Wtime();

    for (int i = 0; i < 10; i++) {
        srand(rank);
        in = 0;
        s = 0;
        for (int i = rank; i < n; i += commsize) {
            double x = getrand() - 1; /* x in [-1, 0] */
            double y = getrand();     /* y in [0, 1] */
            in++;
            s += func(x, y);
        }
        gin = 0;
        MPI_Reduce(&in, &gin, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        gsum = 0.0;
        MPI_Reduce(&s, &gsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    time = MPI_Wtime();
    double tmax = 0;
    MPI_Reduce(&time, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        file = fopen("monte_result.txt", "a");
        fprintf(file, "%f\t%d\n", tmax, commsize);
        fclose(file);
    }
    MPI_Finalize();
    return 0;
}