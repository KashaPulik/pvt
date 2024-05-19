#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <sys/time.h>
#include <utility>

#define N_LOGICAL_CORES 8

int threshold = 1'000;
size_t size = 1'000'000'000;

double wtime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec * 1E-6;
}

void partition(int* v, int& i, int& j, int low, int high)
{
    i = low;
    j = high;
    int pivot = v[(low + high) / 2];
    do {
        while (v[i] < pivot)
            i++;
        while (v[j] > pivot)
            j--;
        if (i <= j) {
            std::swap(v[i], v[j]);
            i++;
            j--;
        }
    } while (i <= j);
}
void quicksort(int* v, int low, int high)
{
    int i, j;
    partition(v, i, j, low, high);
    if (low < j)
        quicksort(v, low, j);
    if (i < high)
        quicksort(v, i, high);
}

void quicksort_tasks(int* v, int low, int high)
{
    int i, j;
    partition(v, i, j, low, high);
    if (high - low < threshold
        || (j - low < threshold || high - i < threshold)) {
        if (low < j)
            quicksort_tasks(v, low, j);
        if (i < high)
            quicksort_tasks(v, i, high);
    } else {
#pragma omp task untied
        {
            quicksort_tasks(v, low, j);
        }
        quicksort_tasks(v, i, high);
    }
}

double run_serial(int* array, int n)
{
    double time = -wtime();
    quicksort(array, 0, n - 1);
    time += wtime();
    return time;
}

double run_parallel(int* array, int n, int n_threads)
{
    double time = -wtime();
#pragma omp parallel num_threads(n_threads)
    {
#pragma omp single
        quicksort_tasks(array, 0, n - 1);
    }
    time += wtime();
    return time;
}

void mix_array(int* array, int size)
{
    srand(0);
    for (int i = 0; i < size; i++)
        std::swap(array[rand() % size], array[rand() % size]);
}

int main()
{
    std::ostringstream filename;
    std::ofstream file;
    int* array;
    double t[N_LOGICAL_CORES + 1];
    for (size = 250'000; size <= 1'000'000; size += 250'000) {
        for (threshold = 200; threshold <= 1200; threshold += 200) {
            filename << size << "and" << threshold << ".dat";
            file.open(filename.str());
            array = new int[size];
            for (size_t i = 0; i < size; i++)
                array[i] = i;
            mix_array(array, size);
            t[1] = run_serial(array, size);
            for (int i = 2; i <= N_LOGICAL_CORES; i += 2) {
                mix_array(array, size);
                t[i] = run_parallel(array, size, i);
            }
            for (int i = 2; i <= N_LOGICAL_CORES; i += 2)
                file << i << "  " << t[1] / t[i] << "\n";
            filename.str("");
            filename.clear();
            delete[] array;
            file.close();
        }
    }
}
