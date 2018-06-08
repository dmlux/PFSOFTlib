
//
//  benchmark_dsoft_runtime_accuracy.cpp
//  PFSOFTlib
//
//   Created by Denis-Michael Lux on 07. Juni 2015.
//
//   This file is part of PFSOFTlib.
//
//   PFSOFTlib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   PFSOFTlib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with PFSOFTlib.  If not, see <http://www.gnu.org/licenses/>.
//

#include <pfsoft>
#include <stdio.h>

using namespace pfsoft;
using namespace FourierTransforms;

void print_error(const DSOFTFourierCoefficients& c1, const DSOFTFourierCoefficients& c2)
{
    // getting errors
    double max_abs_error = 0.0;
    double max_rel_error = 0.0;
    
    // getting BW
    int bandwidth = c1.bandwidth;
    
    // iterate over coefficients
    for (int m = 0; m < bandwidth; ++m)
    {
        for (int n = -m; n <= m; ++n)
        {
            for (int k = -m; k <= m; ++k)
            {
                double abs_error = (c1(m,n,k).re - c2(m,n,k).re) * (c1(m,n,k).re - c2(m,n,k).re) + (c1(m,n,k).im - c2(m,n,k).im) * (c1(m,n,k).im - c2(m,n,k).im);
                double rel_error = ((c1(m,n,k).re - c2(m,n,k).re) * (c1(m,n,k).re - c2(m,n,k).re) + (c1(m,n,k).im - c2(m,n,k).im) * (c1(m,n,k).im - c2(m,n,k).im)) / (c1(m,n,k).re * c1(m,n,k).re + c1(m,n,k).im * c1(m,n,k).im);
                
                if ( abs_error > max_abs_error ) max_abs_error = abs_error;
                if ( rel_error > max_rel_error ) max_rel_error = rel_error;
            }
        }
    }
    
    // normalize errors
    max_abs_error = sqrt (max_abs_error);
    max_rel_error = sqrt (max_rel_error);
    
    // print stuff
    printf("| max abs error:     %.2e\n", max_abs_error);
    printf("| max rel error:     %.2e\n", max_rel_error);
}

int main(int argc, const char** argv)
{
    
    // run benchmark if multithreading is enabled
#ifdef _OPENMP
    int start_threads = 1;
    int BW = 16;
    
    if (argc < 2)
    {
        printf("usage: ./benchmark_sofft_inv_speedup <BANDWIDTH> [<MINIMAL NUMBER OF THREADS>]\n");
        return 1;
    }
    
    if (argc == 2)
    {
        BW = atoi(argv[1]);
    } else if (argc == 3) {
        start_threads = atoi(argv[2]) < 2 ? 2 : atoi(argv[2]);
    }
    
    // To make things fair, we run omp once for startup. This will avoid
    // initialization time later on:
    int max_procs = omp_get_num_procs();
    #pragma omp parallel for num_threads(max_procs)
    for (int i = 0; i < max_procs; i++);
    
    // storage for average runtimes and given number of threads
    double runtimes[omp_get_max_threads() - 1];
    
    // Print to console
    printf("+--------------------------------------------------------------------------------------+\n");
    printf("|                            DSOFT INVERSE SPEEDUP BENCHMARK                           |\n");
    printf("+--------------------------------------------------------------------------------------+\n");
    printf("| FOR BANDWIDTH %i\n", BW);
    printf("| PARALLELIZED WITH %d THREADS\n", omp_get_max_threads());
    printf("+--------------------------------------------------------------------------------------+\n");
    
    // create a grid to fill with values
    grid3D< complex< double > > sample(2 * BW);
    
    // create indices
    unsigned int i, threads;
    
    // creating fourier coefficients container
    DSOFTFourierCoefficients coef(BW);
    DSOFTFourierCoefficients rec_coef(BW);
    
    // generate random coefficients between -1 and 1
    uniform_real_distribution< double > ctx;
    ctx.engine = random_engine::MERSENNE_TWISTER64;
    ctx.min = -1;
    ctx.max = +1;
    
    rand(coef, ctx);
    
    // create sample
    // get reference value of serial implementation
    printf("| Threads:           %d\n", 1);
    
    // compute grid from coefficients
    stopwatch sw_inv = stopwatch::tic();
    IDSOFT(coef, sample, 1);
    double serial_inv_ref = sw_inv.toc();
    printf("| IDSOFT:            %.6fs\n", serial_inv_ref);
    
    // compute grid back to coefficients
    stopwatch sw_for = stopwatch::tic();
    DSOFT(sample, rec_coef, 1);
    double serial_for_ref = sw_for.toc();
    printf("| DSOFT:             %.6fs\n", serial_for_ref);
    
    print_error(coef, rec_coef);
    printf("+--------------------------------------------------------------------------------------+\n");
    
    // run loop run for all number of available threads
    for (threads = start_threads; threads <= omp_get_max_threads(); ++threads)
    {
        printf("| Threads:           %d\n", threads);
        rand(coef, ctx);
        
        // inverse transformation
        sw_inv = stopwatch::tic();
        IDSOFT(coef, sample, threads);
        double inv_runtime = sw_inv.toc();
        printf("| IDSOFT:            %.6fs\n", inv_runtime);
        printf("| Speedup IDSOFT:    %.6fs\n", serial_inv_ref / inv_runtime);
        printf("| Efficiency IDSOFT: %.6fs\n", (serial_inv_ref / inv_runtime) / threads);
        
        // forward transformation
        sw_for = stopwatch::tic();
        DSOFT(sample, rec_coef, threads);
        double for_runtime = sw_for.toc();
        printf("| DSOFT:             %.6fs\n", for_runtime);
        printf("| Speedup DSOFT:     %.6fs\n", serial_for_ref / for_runtime);
        printf("| Efficiency DSOFT:  %.6fs\n", (serial_for_ref / for_runtime) / threads);
        
        print_error(coef, rec_coef);
        printf("+--------------------------------------------------------------------------------------+\n");
    }
    
    return 0;
#else
    printf("Multithreading disabled. Use another compiler that supports OpenMP to run this benchmark!");
    return 1;
#endif
}
