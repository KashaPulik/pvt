#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N_LOGICAL_CORES 8
#define EPS 1e-5
#define N 10000000

const double a = 1, b = 1.2;

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

double getrand(unsigned int* seed)
{
    return (double)rand_r(seed) / RAND_MAX;
}

double func(double x, double y)
{
    return exp(x - y);
}

void integral()
{
    int in = 0;
    unsigned int seed = 0;
    double s = 0;
    for (int i = 0; i < N; i++) {
        double x = -getrand(&seed); /* x in [-1, 0] */
        double y = getrand(&seed);  /* y in [0, 1] */
        in++;
        s += func(x, y);
    }
    double v = in / N;
    double res = v * s / in;
#if 0
    printf("Result: %.12f, n %d\n", res, n);
#endif
}

void integral_omp(int n_threads)
{
    int in = 0;
    double s = 0;
#pragma omp parallel num_threads(n_threads)
    {
        double s_loc = 0;
        int in_loc = 0;
        unsigned int seed = omp_get_thread_num();
#pragma omp for nowait
        for (int i = 0; i < N; i++) {
            double x = -getrand(&seed); /* x in [-1, 0] */
            double y = getrand(&seed);  /* y in [0, 1] */
            in_loc++;
            s_loc += func(x, y);
        }
#pragma omp atomic
        s += s_loc;
#pragma omp atomic
        in += in_loc;
    }
    double v = in / N;
    double res = v * s / in;
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