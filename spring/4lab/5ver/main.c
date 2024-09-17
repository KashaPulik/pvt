#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N_LOGICAL_CORES 8

struct particle {
    float x, y, z;
};

const float G = 6.67e-11;

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

void calculate_forces(struct particle* p, struct particle* f, float* m, int n)
{
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            // Вычисление силы, действующей на тело i со стороны j
            float dist
                    = sqrtf(powf(p[i].x - p[j].x, 2) + powf(p[i].y - p[j].y, 2)
                            + powf(p[i].z - p[j].z, 2));
            float mag = (G * m[i] * m[j]) / powf(dist, 2);
            struct particle dir
                    = {.x = p[j].x - p[i].x,
                       .y = p[j].y - p[i].y,
                       .z = p[j].z - p[i].z};
            // Сумма сил, действующих на тело i
            f[i].x += mag * dir.x / dist;
            f[i].y += mag * dir.y / dist;
            f[i].z += mag * dir.z / dist;
            // Сумма сил, действующих на тело j (симметричность)
            f[j].x -= mag * dir.x / dist;
            f[j].y -= mag * dir.y / dist;
            f[j].z -= mag * dir.z / dist;
        }
    }
}

void move_particles(
        struct particle* p,
        struct particle* f,
        struct particle* v,
        float* m,
        int n,
        float dt)
{
    for (int i = 0; i < n; i++) {
        struct particle dv = {
                .x = f[i].x / m[i] * dt,
                .y = f[i].y / m[i] * dt,
                .z = f[i].z / m[i] * dt,
        };
        struct particle dp = {
                .x = (v[i].x + dv.x / 2) * dt,
                .y = (v[i].y + dv.y / 2) * dt,
                .z = (v[i].z + dv.z / 2) * dt,
        };
        v[i].x += dv.x;
        v[i].y += dv.y;
        v[i].z += dv.z;
        p[i].x += dp.x;
        p[i].y += dp.y;
        p[i].z += dp.z;
        f[i].x = f[i].y = f[i].z = 0;
    }
}

double nbody_serial(int n)
{
    struct particle* p = malloc(sizeof(*p) * n); // Положение частиц (x, y, z)
    struct particle* f = malloc(
            sizeof(*f) * n); // Сила, действующая на каждую частицу (x, y, z)
    struct particle* v = malloc(sizeof(*v) * n); // Скорость частицы (x, y, z)
    float* m = malloc(sizeof(*m) * n); // Масса частицы
    for (int i = 0; i < n; i++) {
        p[i].x = rand() / (float)RAND_MAX - 0.5;
        p[i].y = rand() / (float)RAND_MAX - 0.5;
        p[i].z = rand() / (float)RAND_MAX - 0.5;
        v[i].x = rand() / (float)RAND_MAX - 0.5;
        v[i].y = rand() / (float)RAND_MAX - 0.5;
        v[i].z = rand() / (float)RAND_MAX - 0.5;
        m[i] = rand() / (float)RAND_MAX * 10 + 0.01;
        f[i].x = f[i].y = f[i].z = 0;
    }
    double dt = 1e-5;
    double time = -wtime();
    for (double t = 0; t <= 1; t += dt) { // Цикл по времени (модельному)
        calculate_forces(p, f, m, n); // Вычисление сил – O(N^2)
        move_particles(p, f, v, m, n, dt); // Перемещение тел O(N)
    }
    time += wtime();
    free(m);
    free(v);
    free(f);
    free(p);
    return time;
}

void calculate_forces_parallel(
        struct particle* p, struct particle* f[], float* m, int n)
{
    int tid = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    for (int i = 0; i < n; i++) {
        f[tid][i].x = 0;
        f[tid][i].y = 0;
        f[tid][i].z = 0;
    }
#pragma omp for schedule(dynamic, 8)
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            float dist
                    = sqrtf(powf(p[i].x - p[j].x, 2) + powf(p[i].y - p[j].y, 2)
                            + powf(p[i].z - p[j].z, 2));
            float mag = (G * m[i] * m[j]) / powf(dist, 2);
            struct particle dir
                    = {.x = p[j].x - p[i].x,
                       .y = p[j].y - p[i].y,
                       .z = p[j].z - p[i].z};
            f[tid][i].x += mag * dir.x / dist;
            f[tid][i].y += mag * dir.y / dist;
            f[tid][i].z += mag * dir.z / dist;
            f[tid][j].x -= mag * dir.x / dist;
            f[tid][j].y -= mag * dir.y / dist;
            f[tid][j].z -= mag * dir.z / dist;
        }
    }
#pragma omp single // Итоговый вектор сил сформируем в первой строке – f[0][i]
    {
        for (int i = 0; i < n; i++) {
            for (int tid = 1; tid < nthreads; tid++) {
                f[0][i].x += f[tid][i].x;
                f[0][i].y += f[tid][i].y;
                f[0][i].z += f[tid][i].z;
            }
        }
    }
}

void move_particles_parallel(
        struct particle* p,
        struct particle* f[],
        struct particle* v,
        float* m,
        int n,
        float dt)
{
#pragma omp for
    for (int i = 0; i < n; i++) {
        struct particle dv = {
                .x = f[0][i].x / m[i] * dt,
                .y = f[0][i].y / m[i] * dt,
                .z = f[0][i].z / m[i] * dt,
        };
        struct particle dp = {
                .x = (v[i].x + dv.x / 2) * dt,
                .y = (v[i].y + dv.y / 2) * dt,
                .z = (v[i].z + dv.z / 2) * dt,
        };
        v[i].x += dv.x;
        v[i].y += dv.y;
        v[i].z += dv.z;
        p[i].x += dp.x;
        p[i].y += dp.y;
        p[i].z += dp.z;
        // f[i].x = f[i].y = f[i].z = 0;
    }
}

double nbody_parallel(int n, int n_threads)
{
    struct particle* p = malloc(sizeof(*p) * n); // Положение частиц (x, y, z)
    struct particle* f[n_threads];
    for (int i = 0; i < n_threads; i++) {
        f[i] = malloc(sizeof(struct particle) * n);
        for (int j = 0; j < n; j++) {
            f[i][j].x = 0;
            f[i][j].y = 0;
            f[i][j].z = 0;
        }
    }
    struct particle* v = malloc(sizeof(*v) * n); // Скорость частицы (x, y, z)
    float* m = malloc(sizeof(*m) * n); // Масса частицы
    for (int i = 0; i < n; i++) {
        p[i].x = rand() / (float)RAND_MAX - 0.5;
        p[i].y = rand() / (float)RAND_MAX - 0.5;
        p[i].z = rand() / (float)RAND_MAX - 0.5;
        v[i].x = rand() / (float)RAND_MAX - 0.5;
        v[i].y = rand() / (float)RAND_MAX - 0.5;
        v[i].z = rand() / (float)RAND_MAX - 0.5;
        m[i] = rand() / (float)RAND_MAX * 10 + 0.01;
    }
    double dt = 1e-5;
    double time = -wtime();
#pragma omp parallel num_threads(n_threads)
    {
        for (double t = 0; t <= 1; t += dt) { // Цикл по времени (модельному)
            calculate_forces_parallel(p, f, m, n); // Вычисление сил – O(N^2)
            move_particles_parallel(p, f, v, m, n, dt); // Перемещение тел O(N)
        }
    }
    time += wtime();
    free(m);
    free(v);
    free(p);
    for (int i = 0; i < n_threads; i++)
        free(f[i]);
    return time;
}

int main(int argc, char* argv[])
{
    int n = (argc > 1) ? atoi(argv[1]) : 10;
    double t[N_LOGICAL_CORES + 1];
    FILE* file;

    t[1] = nbody_serial(n);
    printf("%d threads: %f seconds\n", 1, t[1]);
    for (int i = 2; i <= N_LOGICAL_CORES; i += 2) {
        t[i] = nbody_parallel(n, i);
        printf("%d threads: %f seconds\n", i, t[i]);
    }
    file = fopen("data.dat", "w");
    for (int i = 2; i <= N_LOGICAL_CORES; i += 2) {
        fprintf(file, "%d    %f\n", i, t[1] / t[i]);
    }
    fclose(file);
    return 0;
}