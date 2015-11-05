//
//  benchmark_sofft_inv_speedup.cpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 22.07.15.
//
//

#include <uzlmath>
#include <stdio.h>

using namespace uzlmath;
using namespace FourierTransforms;

int main(int argc, const char** argv)
{
    
    // run benchmark if multithreading is enabled
#ifdef _OPENMP
    if (argc < 4)
    {
        printf("usage: ./benchmark_sofft_inv_speedup <MIN BANDWIDTH> <MAX BANDWIDTH> <RUNS PER BANDWDITH>\n");
        return 1;
    }
    
    int START_BW = atoi(argv[1]);
    int MAX_BW = atoi(argv[2]);
    int LOOP_R = atoi(argv[3]);
    
    // To make things fair, we run omp once for startup. This will avoid
    // initialization time later on:
    int max_procs = omp_get_num_procs();
    #pragma omp parallel for num_threads(max_procs)
    for (int i = 0; i < max_procs; i++);
    
    // storage for average runtimes and given number of threads
    double runtimes[omp_get_max_threads() - 1];
    
    // write to file
    FILE* fp  = fopen("benchmark_DSOFT_inv_speedup.txt", "w");
    FILE* fp2 = fopen("DSOFT_inverse_speedup.dat", "w");
    FILE* fp3 = fopen("DSOFT_inverse_speedup_runtimes.dat", "w");
    
    // Print to console
    printf("+--------------------------------------------------------------------------------------+\n");
    printf("|                            DSOFT INVERSE SPEEDUP BENCHMARK                           |\n");
    printf("+--------------------------------------------------------------------------------------+\n");
    printf("| FROM BANDWIDTH %i TO %i WITH %i LOOP RUNS PER BANDWIDTH\n", START_BW, MAX_BW, LOOP_R);
    printf("| PARALLELIZED WITH %d THREADS\n", omp_get_max_threads());
    
    printf("+=====+============+");
    for (int i = 0; i < omp_get_max_threads() - 1; ++i) { printf("==================+"); }
    printf("\n");
    
    printf("|  B  | t (serial) |");
    for (int i = 0; i < omp_get_max_threads() - 1; ++i) { printf(" speedup %d cores  |", i + 2); }
    printf("\n");
    
    printf("+=====+============+");
    for (int i = 0; i < omp_get_max_threads() - 1; ++i) { printf("==================+"); }
    printf("\n");
    
    // Print to file
    fprintf(fp, "+--------------------------------------------------------------------------------------+\n");
    fprintf(fp, "|                            DSOFT INVERSE SPEEDUP BENCHMARK                           |\n");
    fprintf(fp, "+--------------------------------------------------------------------------------------+\n");
    fprintf(fp, "| FROM BANDWIDTH %i TO %i WITH %i LOOP RUNS PER BANDWIDTH\n", START_BW, MAX_BW, LOOP_R);
    fprintf(fp, "| PARALLELIZED WITH %d THREADS\n", omp_get_max_threads());
    
    fprintf(fp, "+=====+============+");
    for (int i = 0; i < omp_get_max_threads() - 1; ++i) { fprintf(fp, "==================+"); }
    fprintf(fp, "\n");
    
    fprintf(fp, "|  B  | t (serial) |");
    for (int i = 0; i < omp_get_max_threads() - 1; ++i) { fprintf(fp, " speedup %d cores  |", i + 2); }
    fprintf(fp, "\n");
    
    fprintf(fp, "+=====+============+");
    for (int i = 0; i < omp_get_max_threads() - 1; ++i) { fprintf(fp, "==================+"); }
    fprintf(fp, "\n");
    
    // print labels to file
    fprintf(fp2, "bandwidth\tserial\t");
    for (int i = 0; i < omp_get_max_threads() - 1; ++i) { fprintf(fp2, "c%d\ttc%d\t", i + 2, i + 2); }
    fprintf(fp2, "\n");
    
    fprintf(fp3, "bandwidth\tserial\tthreads\t");
    for (int i = 0; i < LOOP_R; ++i)
    {
        fprintf(fp3, "c%i\t", i+1);
    }
    fprintf(fp3, "\n");
    
    // run benchmark up to MAX_BW times
    for (int bandwidth = START_BW; bandwidth <= MAX_BW; ++bandwidth)
    {
        
        // create a grid to fill with values
        grid3D< complex< double > > sample(2 * bandwidth);
        
        // create indices
        unsigned int i, threads;
        
        // creating fourier coefficients container
        DSOFTFourierCoefficients coef(bandwidth);
        
        // generate random coefficients between -1 and 1
        uniform_real_distribution< double > ctx;
        ctx.engine = random_engine::MERSENNE_TWISTER64;
        ctx.min = -1;
        ctx.max = +1;
        
        rand(coef, ctx);
        
        // create sample
        // get reference value of serial implementation
        double serial_ref = 0;
        for (int i = 0; i < LOOP_R; ++i)
        {
            stopwatch sw = stopwatch::tic();
            IDSOFT(coef, sample, 1);
            serial_ref += sw.toc();
        }
        serial_ref /= LOOP_R;
        
        // reset times
        for (i = 0; i < omp_get_max_threads(); ++i)
        {
            runtimes[i] = 0;
        }
        
        // run loop run for all number of available threads
        for (threads = 2; threads <= omp_get_max_threads(); ++threads)
        {
            fprintf(fp3, "%i\t%3.6f\t%i\t", bandwidth, serial_ref, threads);
            // run the needed amount of loopruns
            for (i = 0; i < LOOP_R; ++i)
            {
                // perform forward DSOFT transform
                // and stop time
                stopwatch sw = stopwatch::tic();
                IDSOFT(coef, sample, threads);
                double time  = sw.toc();
                
                fprintf(fp3, "%3.6f\t", time);
                
                runtimes[threads - 2] += time;
            }
            fprintf(fp3, "\n");
        }
        
        fprintf(fp2, "%3d\t%2.6f\t", bandwidth, serial_ref);
        for (int i = 0; i < omp_get_max_threads() - 1; ++i)
        {
            fprintf(fp2, "%2.2f\t%2.6f\t", (serial_ref / (runtimes[i] / LOOP_R)), (runtimes[i] / LOOP_R));
        }
        fprintf(fp2, "\n");
        
        // print information
        printf("| %3d | %2.6fs  | ", bandwidth, serial_ref);
        for (i = 0; i < omp_get_max_threads() - 1; ++i)
        {
            printf("%2.2f (%2.6fs) | ", (serial_ref / (runtimes[i] / LOOP_R)), (runtimes[i] / LOOP_R));
        }
        printf("\n");
        
        // print info to file
        fprintf(fp, "| %3d | %2.6fs  | ", bandwidth, serial_ref);
        for (i = 0; i < omp_get_max_threads() - 1; ++i)
        {
            fprintf(fp, "%2.2f (%2.6fs) | ", (serial_ref / (runtimes[i] / LOOP_R)), (runtimes[i] / LOOP_R));
        }
        fprintf(fp, "\n");
    }
    
    printf("+=====+============+");
    for (int i = 0; i < omp_get_max_threads() - 1; ++i) { printf("==================+"); }
    printf("\n");
    
    fprintf(fp, "+=====+============+");
    for (int i = 0; i < omp_get_max_threads() - 1; ++i) { fprintf(fp, "==================+"); }
    fprintf(fp, "\n");
    
    // close files
    fclose( fp  );
    fclose( fp2 );
    fclose( fp3 );
    
    return 0;
#else
    printf("Multithreading disabled. Use another compiler that supports OpenMP to run this benchmark!");
    return 1;
#endif
}