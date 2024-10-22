#include <math.h>
#include <mpi.h>
#include <stdio.h>

const double PI = 3.14159265358979323846;

double func(double x)
{
    return sqrt(x * (3 - x)) / (x + 1);
}

const double a = 1;
const double b = 1.2;
const double eps = 1e-14;
const int n0 = 100;

int main(int argc, char** argv)
{
    FILE* file;

    int commsize, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double time = MPI_Wtime();
    int n = n0, k;
    double sq[2], delta = 1;
    for (int i = 0; i < 10; i++) {
        delta = 1;
        for (k = 0; delta > eps; n *= 2, k ^= 1) {
            int points_per_proc = n / commsize;
            int lb = rank * points_per_proc;
            int ub = (rank == commsize - 1) ? (n - 1) : (lb + points_per_proc - 1);
            double h = (b - a) / n;
            double s = 0.0;
            for (int i = lb; i <= ub; i++)
                s += func(a + h * (i + 0.5));
            MPI_Allreduce(&s, &sq[k], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            sq[k] *= h;
            if (n > n0)
                delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
        }
    }
    time = MPI_Wtime() - time;
    double tmax = 0;
    MPI_Reduce(&time, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        file = fopen("midpoint_result.txt", "a");
        fprintf(file, "%f\t%d\n", tmax, commsize);
        fclose(file);
    }
    MPI_Finalize();
    return 0;
}
