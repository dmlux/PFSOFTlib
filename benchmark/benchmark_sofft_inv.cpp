//
//  benchmark_sofft_inv.cpp
//  PDSOFTlib
//
//   Created by Denis-Michael Lux on 05. November 2015.
//
//   This file is part of PDSOFTlib.
//
//   PDSOFTlib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   PDSOFTlib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with PDSOFTlib.  If not, see <http://www.gnu.org/licenses/>.
//

#include <pdsoft>
#include <stdio.h>

using namespace pdsoft;
using namespace FourierTransforms;

// Main method
int main(int argc, const char** argv)
{
    if (argc < 4)
    {
        printf("usage: ./benchmark_sofft_inv <MIN BANDWIDTH> <MAX BANDWIDTH> <RUNS PER BANDWDITH>\n");
        return 1;
    }
    
    int START_BW = atoi(argv[1]);
    int MAX_BW = atoi(argv[2]);
    int LOOP_R = atoi(argv[3]);
    
    // To make things fair, we run omp once for startup. This will avoid
    // initialization time later on:
#ifdef _OPENMP
    int max_procs = omp_get_num_procs();
    #pragma omp parallel for num_threads(max_procs)
    for (int i = 0; i < max_procs; i++);
#endif
    
    // write to file
    FILE* fp  = fopen("benchmark_DSOFT_inv.txt", "w");
    FILE* fp2 = fopen("DSOFT_inverse.dat", "w");

    // write output to file "benchmark_sofft_inv.txt"
#ifdef _OPENMP
    printf(     "+-----------------------------------------------------------------------------------------------------------------------------+\n");
    printf(     "|                                                    DSOFT INVERSE BENCHMARK                                                  |\n");
    printf(     "+-----------------------------------------------------------------------------------------------------------------------------+\n");
    printf(     "| FROM BANDWIDTH %i TO %i WITH %i LOOP RUNS PER BANDWIDTH\n", START_BW, MAX_BW, LOOP_R);
    printf(     "| PARALLELIZED WITH %d THREADS\n", omp_get_max_threads());
    printf(     "+=====+===========+===================================+===================================+===========+==========+============+\n");
    printf(     "|  B  | average   | fastest run (dif. to avg / %%dif)  | slowest run (dif. to avg / %%dif)  | serial    | speedup  | efficiency |\n");
    printf(     "+=====+===========+===================================+===================================+===========+==========+============+\n");
    
    fprintf(fp, "+-----------------------------------------------------------------------------------------------------------------------------+\n");
    fprintf(fp, "|                                                    DSOFT INVERSE BENCHMARK                                                  |\n");
    fprintf(fp, "+-----------------------------------------------------------------------------------------------------------------------------+\n");
    fprintf(fp, "| FROM BANDWIDTH %i TO %i WITH %i LOOP RUNS PER BANDWIDTH\n", START_BW, MAX_BW, LOOP_R);
    fprintf(fp, "| PARALLELIZED WITH %d THREADS\n", omp_get_max_threads());
    fprintf(fp, "+=====+===========+===================================+===================================+===========+==========+============+\n");
    fprintf(fp, "|  B  | average   | fastest run (dif. to avg / %%dif)  | slowest run (dif. to avg / %%dif)  | serial    | speedup  | efficiency |\n");
    fprintf(fp, "+=====+===========+===================================+===================================+===========+==========+============+\n");
    
    fprintf(fp2, "bandwidth\truntime\tserial\tspeedup\tefficiency\n");
#else
    printf(     "+-----------------------------------------------------------------------------------------+\n");
    printf(     "|                                  DSOFT INVERSE BENCHMARK                                |\n");
    printf(     "+-----------------------------------------------------------------------------------------+\n");
    printf(     "| FROM BANDWIDTH %i TO %i WITH %i LOOP RUNS PER BANDWIDTH\n", START_BW, MAX_BW, LOOP_R);
    printf(     "+=====+===========+===================================+===================================+\n");
    printf(     "|  B  | average   | fastest run (dif. to avg / %%dif)  | slowest run (dif. to avg / %%dif)  |\n");
    printf(     "+=====+===========+===================================+===================================+\n");
    
    fprintf(fp, "+-----------------------------------------------------------------------------------------+\n");
    fprintf(fp, "|                                  DSOFT INVERSE BENCHMARK                                |\n");
    fprintf(fp, "+-----------------------------------------------------------------------------------------+\n");
    fprintf(fp, "| FROM BANDWIDTH %i TO %i WITH %i LOOP RUNS PER BANDWIDTH\n", START_BW, MAX_BW, LOOP_R);
    fprintf(fp, "+=====+===========+===================================+===================================+\n");
    fprintf(fp, "|  B  | average   | fastest run (dif. to avg / %%dif)  | slowest run (dif. to avg / %%dif)  |\n");
    fprintf(fp, "+=====+===========+===================================+===================================+\n");
    
    fprintf(fp2, "bandwidth\truntime\n");
#endif
    
    // loop over all bandwidth up to MAX_BW
    for (unsigned int bandwidth = START_BW; bandwidth <= MAX_BW; ++bandwidth)
    {
        
        // create a grid to fill with values
        grid3D< complex< double > > sample(2 * bandwidth);
        
        // create indices
        unsigned int i;
        
        // creating fourier coefficients container
        DSOFTFourierCoefficients coef(bandwidth);
        
        // generate random coefficients between -1 and 1
        uniform_real_distribution< double > ctx;
        ctx.engine = random_engine::MERSENNE_TWISTER64;
        ctx.min = -1;
        ctx.max = +1;
        
        rand(coef, ctx);
        
        // min and max exec tiems
        double min, max;
        
        // variable for execution times
        double times = 0;
        
#ifdef _OPENMP
        // get reference value of serial implementation
        double serial_ref = 0;
        for (int i = 0; i < LOOP_R; ++i)
        {
            stopwatch sw = stopwatch::tic();
            IDSOFT(coef, sample, 1);  // setting threads explicitly to 1
            serial_ref += sw.toc();
        }
        serial_ref /= LOOP_R;
#endif
        
        for (i = 0; i < LOOP_R; ++i)
        {
            // perform inverse DSOFT transform
            // and stop time
            stopwatch sw = stopwatch::tic();
            IDSOFT(coef, sample);
            double time  = sw.toc();
            
            // add to sum of time for current bandwidth
            times += time;
            
            // if first run then set min and max to the
            // first runtime
            if (i == 0)
            {
                min = time;
                max = time;
            }
            // store the current value if it is smaller than
            // the stored minimum (then it must be the new
            // minimum), or if it is bigger than max (then it
            // must be the new maximum).
            else
            {
                min = (time < min ? time : min);
                max = (time > max ? time : max);
            }
        }
        
        // get average execution time
        double avg = times / LOOP_R;
        
        // get fastest and slowest run-ratios
        double min_ratio = (avg - min)/avg * 100;
        double max_ratio = (max - avg)/avg * 100;
        
        // format ratios to string for right alignment
        char min_rat[7], max_rat[7];
        snprintf(min_rat, 7, "%3.3f", min_ratio);
        snprintf(max_rat, 7, "%3.3f", max_ratio);
        
        // print information
        printf("| %3i |", bandwidth);                                                // bandwidth
        printf(" %2.6fs |", avg);                                                    // average runtime for given bandwidth
        printf(" %2.6fs (-%2.6fs / -%6s%%) |", min, (avg - min), min_rat);           // fastest run and its difference to average
        printf(" %2.6fs (+%2.6fs / +%6s%%) |", max, (max - avg), max_rat);           // slowest run and its difference to average
        
        // write to file
        fprintf(fp, "| %3i |", bandwidth);                                           // bandwidth
        fprintf(fp, " %2.6fs |", avg);                                               // average runtime for given bandwidth
        fprintf(fp, " %2.6fs (-%2.6fs / -%6s%%) |", min, (avg - min), min_rat);      // fastest run and its difference to average
        fprintf(fp, " %2.6fs (+%2.6fs / +%6s%%) |", max, (max - avg), max_rat);      // slowest run and its difference to average
        
#ifdef _OPENMP
        printf(     " %2.6fs |", serial_ref);
        printf(     "   %2.2f   |", (serial_ref / max));
        printf(     "    %2.2f    |\n", (serial_ref / (omp_get_max_threads() * max)));
        
        fprintf(fp, " %2.6fs |", serial_ref);
        fprintf(fp, "   %2.2f   |", (serial_ref / max));
        fprintf(fp, "    %2.2f    |\n", (serial_ref / (omp_get_max_threads() * max)));
        
        fprintf(fp2, "%d\t\t%15f\t\t%15f\t\t%15f\t\t%15f\n", bandwidth, avg, serial_ref, (serial_ref / max), (serial_ref / (omp_get_max_threads() * max)));
#else
        printf("\n");
        fprintf(fp, "\n");
        
        fprintf(fp2, "%d\t\t%15f\n", bandwidth, avg);
#endif
        
    }
    
#ifdef _OPENMP
    printf(     "+=====+===========+===================================+===================================+===========+==========+============+\n");
    fprintf(fp, "+=====+===========+===================================+===================================+===========+==========+============+\n");
#else
    printf(     "+=====+===========+===================================+===================================+\n");
    fprintf(fp, "+=====+===========+===================================+===================================+\n");
#endif
    
    // close files
    fclose(fp);
    fclose(fp2);
}
