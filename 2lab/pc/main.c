#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N_LOGICAL_CORES 24

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

void matrix_vector_product_omp(
        double* a, double* b, double* c, int m, int n, int n_threads)
{
#pragma omp parallel num_threads(n_threads)
    {
        int nthreads = omp_get_num_threads();
        int threadid = omp_get_thread_num();
        int items_per_thread = m / nthreads;
        int lb = threadid * items_per_thread;
        int ub = (threadid == nthreads - 1) ? (m - 1)
                                            : (lb + items_per_thread - 1);
        for (int i = lb; i <= ub; i++) {
            c[i] = 0.0;
            for (int j = 0; j < n; j++)
                c[i] += a[i * n + j] * b[j];
        }
    }
}

void matrix_vector_product(double* a, double* b, double* c, int m, int n)
{
    for (int i = 0; i < m; i++) {
        c[i] = 0.0;
        for (int j = 0; j < n; j++)
            c[i] += a[i * n + j] * b[j];
    }
}

double run_serial(int m, int n)
{
    double *a, *b, *c;
    a = malloc(sizeof(*a) * m * n);
    b = malloc(sizeof(*b) * n);
    c = malloc(sizeof(*c) * m);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            a[i * n + j] = i + j;
    }
    for (int j = 0; j < n; j++)
        b[j] = j;
    double t = wtime();
    matrix_vector_product(a, b, c, m, n);
    t = wtime() - t;
    printf("%d x %d elements\n", m, n);
    printf("Elapsed time (serial): %.6f sec.\n", t);
    free(a);
    free(b);
    free(c);
    return t;
}

double run_parallel(int m, int n, int n_threads)
{
    double *a, *b, *c;
    // Allocate memory for 2-d array a[m, n]
    a = malloc(sizeof(*a) * m * n);
    b = malloc(sizeof(*b) * n);
    c = malloc(sizeof(*c) * m);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            a[i * n + j] = i + j;
    }
    for (int j = 0; j < n; j++)
        b[j] = j;
    double t = wtime();
    matrix_vector_product_omp(a, b, c, m, n, n_threads);
    t = wtime() - t;
    printf("Elapsed time (parallel: %d threads): %.6f sec.\n", n_threads, t);
    free(a);
    free(b);
    free(c);
    return t;
}

int main(int argc, char** argv)
{
    int n[] = {15000, 20000, 25000};
    double t[N_LOGICAL_CORES + 1];
    FILE* file;

    t[1] = run_serial(n[0], n[0]);
    for (int i = 2; i <= N_LOGICAL_CORES; i += 2)
        t[i] = run_parallel(n[0], n[0], i);
    file = fopen("15000.dat", "w");
    for (int i = 2; i <= N_LOGICAL_CORES; i += 2) {
        fprintf(file, "%d    %f\n", i, t[1] / t[i]);
    }
    fclose(file);

    t[1] = run_serial(n[1], n[1]);
    for (int i = 2; i <= N_LOGICAL_CORES; i += 2)
        t[i] = run_parallel(n[1], n[1], i);
    file = fopen("20000.dat", "w");
    for (int i = 2; i <= N_LOGICAL_CORES; i += 2) {
        fprintf(file, "%d    %f\n", i, t[1] / t[i]);
    }
    fclose(file);

    t[1] = run_serial(n[2], n[2]);
    for (int i = 2; i <= N_LOGICAL_CORES; i += 2)
        t[i] = run_parallel(n[2], n[2], i);
    file = fopen("25000.dat", "w");
    for (int i = 2; i <= N_LOGICAL_CORES; i += 2) {
        fprintf(file, "%d    %f\n", i, t[1] / t[i]);
    }
    fclose(file);

    return 0;
}