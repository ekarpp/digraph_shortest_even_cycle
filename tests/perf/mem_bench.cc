/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <omp.h>
#include <stdint.h>

using namespace std;

// 500 mil = 4 GB (500mil * 8B = 4 bil)
#define N (1ull << 29)
#define GIBS (N*8 / (1ull << 30))

int main(void)
{
    double start, end, delta;

    vector<uint64_t> vec(N);
    uint64_t val = time(nullptr);

    start = omp_get_wtime();
    for (uint64_t i = 0; i < N; i++)
        vec[i] = val;
    end = omp_get_wtime();

    delta = end - start;
    cout << "wrote sequentially " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    uint64_t sum = 0;
    start = omp_get_wtime();
    for (uint64_t i = 0; i < N; i++)
        sum += vec[i];
    end = omp_get_wtime();
    cout << sum << endl;

    delta = end - start;
    cout << "read sequentially " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    val = time(nullptr);
    start = omp_get_wtime();
    #pragma omp parallel for
    for (uint64_t i = 0; i < N; i++)
        vec[i] = val;
    end = omp_get_wtime();

    delta = end - start;
    cout << "wrote parallelly " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    uint64_t sums[128];
    sum = 0;
    int cores = omp_get_num_procs();
    uint64_t block_size = N / cores;
    start = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < cores; i++)
    {
        uint64_t csum = 0;
        uint64_t start = block_size * i;
        uint64_t end = block_size * (i + 1) - 1;
        for (uint64_t j = start; j < end; j++)
            csum += vec[j];
        sums[i] = csum;
    }
    for (int i = 0; i < cores; i++)
        sum += sums[i];
    end = omp_get_wtime();
    cout << sum << endl;

    delta = end - start;
    cout << "read parallelly " << GIBS
         << " GiB in " << delta << " seconds, "
         << GIBS / delta << " GiB / s." << endl;

    return 0;
}
