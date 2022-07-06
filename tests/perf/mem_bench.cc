/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#include <iostream>
#include <vector>
#include <omp.h>
#include <stdint.h>

using namespace std;

// 500 mil = 4 GB (500mil * 8B = 4 bil)
#define N 500000000ull
#define GBS (N*8 / 1e9)

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
    cout << "wrote sequentially " << GBS
         << " GB in " << delta << " seconds." << endl;

    uint64_t sum = 0;
    start = omp_get_wtime();
    for (uint64_t i = 0; i < N; i++)
        sum += vec[i];
    end = omp_get_wtime();
    cout << sum << endl;

    delta = end - start;
    cout << "read sequentially " << GBS
         << " GB in " << delta << " seconds." << endl;

    val = time(nullptr);
    start = omp_get_wtime();
    #pragma omp parallel for
    for (uint64_t i = 0; i < N; i++)
        vec[i] = val;
    end = omp_get_wtime();

    delta = end - start;
    cout << "wrote parallelly " << GBS
         << " GB in " << delta << " seconds." << endl;

    uint64_t sums[128];
    sum = 0;
    int threads = omp_get_num_threads();
    uint64_t block_size = N / threads;
    start = omp_get_wtime();
    for (int i = 0; i < threads; i++)
    {
        uint64_t tsum = 0;
        uint64_t start = block_size * i;
        uint64_t end = block_size * (i + 1) - 1;
        for (uint64_t j = start; j < end; j++)
            tsum += vec[j];
        sums[i] = tsum;
    }
    for (int i = 0; i < threads; i++)
        sum += sums[i];
    end = omp_get_wtime();
    cout << sum << endl;

    delta = end - start;
    cout << "read parallelly " << GBS
         << " GB in " << delta << " seconds." << endl;

    return 0;
}
