#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N_LOGICAL_CORES 8
#define EPS 1e-5
#define N_CUTS 100000000

const double a = 1, b = 1.2;

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

double func(double x)
{
    return sqrt(x * (3.0 - x)) / (x + 1.0);
}

void integral()
{
    size_t n = N_CUTS, k;
    double sq[2], delta = 1;
    for (k = 0; delta > EPS; n *= 2, k ^= 1) {
        double h = (b - a) / n;
        double s = 0.0;
        for (size_t i = 0; i < n; i++)
            s += func(a + h * (i + 0.5));
        sq[k] = s * h;
        if (n > N_CUTS)
            delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
    }
#if 0
    printf("n=%ld i=%ld sq=%.12f delta=%.12f\n", n, k, sq[k], delta);
#endif
}

void integral_omp(int n_threads)
{
    double sq[2];
#pragma omp parallel num_threads(n_threads)
    {
        size_t n = N_CUTS, k;
        double delta = 1;
        for (k = 0; delta > EPS; n *= 2, k ^= 1) {
            double h = (b - a) / n;
            double s = 0.0;
            sq[k] = 0;
// Ждем пока все потоки закончат обнуление sq[k], s
#pragma omp barrier
#pragma omp for nowait
            for (int i = 0; i < n; i++)
                s += func(a + h * (i + 0.5));
#pragma omp atomic
            sq[k] += s * h;
// Ждем пока все потоки обновят sq[k]
#pragma omp barrier
            if (n > N_CUTS)
                delta = fabs(sq[k] - sq[k ^ 1]) / 3.0;
#if 0
printf("n=%ld i=%ld sq=%.12f delta=%.12f\n", n, k, sq[k], delta);
#endif
        }
    }
}

double run_serial()
{
    double t = -wtime();
    integral();
    t += wtime();
    printf("Elapsed time (serial): %.6f sec.\n", t);
    return t;
}

double run_parallel(int n_threads)
{
    double t = -wtime();
    integral_omp(n_threads);
    t += wtime();
    printf("Elapsed time (parallel: %d threads): %.6f sec.\n", n_threads, t);
    return t;
}

int main(int argc, char** argv)
{
    double t[N_LOGICAL_CORES + 1];
    FILE* file;

    t[1] = run_serial();
    for (int i = 2; i <= N_LOGICAL_CORES; i++)
        t[i] = run_parallel(i);
    file = fopen("data.dat", "w");
    for (int i = 1; i <= N_LOGICAL_CORES; i++) {
        fprintf(file, "%d    %f\n", i, t[1] / t[i]);
    }
    fclose(file);
    return 0;
}