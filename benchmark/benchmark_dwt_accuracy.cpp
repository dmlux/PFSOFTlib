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

int main(int argc, const char** argv)
{
    if (argc < 3)
    {
        printf("usage: ./benchmark_dwt_for_accuracy <B> <RUNS>\n");
        return 1;
    }
    
    int B    = atoi(*(argv + 1));
    int runs = atoi(*(argv + 2));
    
    printf("+-----------------------------------------------------------------------+\n");
    printf("|                         BENCHMARK DWT ACCURACY                        |\n");
    printf("+-----------------------------------------------------------------------+\n");
    
    printf("+------+--------------+--------------+--------------+-----------------+\n");
    printf("|  BW  | M=0, M'=0    | M=BW/2, M'=0 | M=M'=BW/2    | %d iterations |\n", runs);
    printf("+------+--------------+--------------+--------------+-----------------+\n");
    
    // run all power of 2 bandwidths
    for (int bw = 2; bw <= B; bw *= 2)
    {
        
        // create weights for the given bandwidth
        vector< double > weights(2 * bw);
        DWT::quadrature_weights(weights);
    
        // create weigthed wigner matrix and wigner matrix
        matrix< double > dw(bw, 2 * bw);
        DWT::weighted_wigner_d_matrix(dw, bw, 0, 0, weights);
        
        matrix< double > dt(bw, 2 * bw);
        DWT::wigner_d_matrix(dw, bw, 0, 0);
        dt.transpose();
        
        // relative and absolute errors
        double relative[3];
        double absolute[3];
        
        // init errors
        for (int i = 0; i < 3; ++i)
        {
            relative[i] = 0;
            absolute[i] = 0;
        }
        
        // run several runs and take average values
        for (int i = 0; i < runs; ++i)
        {
            
            // create random coefficients
            vector< complex< double > > fh(dt.cols, vector< complex< double > >::type::COLUMN);
            rand(fh, -1, 1);
            
            // inverse DWT
            vector< complex< double > > s = dt * fh;
            
            // forward DWT
            vector< complex< double > > gh = dw * s;
            
            // get difference vector
            vector< complex< double > > dif = gh - fh;
            
            // get length of error vector
            double abs = 0;
            double org = 0;
            
            for (int j = 0; j < dif.size; ++j)
            {
                abs += dif[j].abs() * dif[j].abs();
                org += fh[j].abs() * fh[j].abs();
            }
            
            relative[0] += sqrt(abs) / sqrt(org);
            absolute[0] += sqrt(abs);
        }
        
        dw = matrix< double >(bw/2, 2 * bw);
        DWT::weighted_wigner_d_matrix(dw, bw, bw/2, 0, weights);
        
        dt = matrix< double >(bw/2, 2 * bw);
        DWT::wigner_d_matrix(dt, bw, bw/2, 0);
        dt.transpose();
        
        // run several runs and take average values
        for (int i = 0; i < runs; ++i)
        {
            
            // create random coefficients
            vector< complex< double > > fh(dt.cols, vector< complex< double > >::type::COLUMN);
            rand(fh, -1, 1);
            
            // inverse DWT
            vector< complex< double > > s = dt * fh;
            
            // forward DWT
            vector< complex< double > > gh = dw * s;
            
            // get difference vector
            vector< complex< double > > dif = gh - fh;
            
            // get length of error vector
            double abs = 0;
            double org = 0;
            for (int j = 0; j < dif.size; ++j)
            {
                abs += dif[j].abs() * dif[j].abs();
                org += fh[j].abs() * fh[j].abs();
            }
            
            relative[1] += sqrt(abs) / sqrt(org);
            absolute[1] += sqrt(abs);
        }
        
        dw = matrix< double >(bw/2, 2 * bw);
        DWT::weighted_wigner_d_matrix(dw, bw, bw/2, bw/2, weights);
        
        dt = matrix< double >(bw/2, 2 * bw);
        DWT::wigner_d_matrix(dt, bw, bw/2, bw/2);
        dt.transpose();
        
        // run several runs and take average values
        for (int i = 0; i < runs; ++i)
        {
            
            // create random coefficients
            vector< complex< double > > fh(dt.cols, vector< complex< double > >::type::COLUMN);
            rand(fh, -1, 1);
            
            // inverse DWT
            vector< complex< double > > s = dt * fh;
            
            // forward DWT
            vector< complex< double > > gh = dw * s;
            
            // get difference vector
            vector< complex< double > > dif = gh - fh;
            
            // get length of error vector
            double abs = 0;
            double org = 0;
            for (int j = 0; j < dif.size; ++j)
            {
                abs += dif[j].abs() * dif[j].abs();
                org += fh[j].abs() * fh[j].abs();
            }
            
            relative[2] += sqrt(abs) / sqrt(org);
            absolute[2] += sqrt(abs);
        }
        
        // get average value
        for (int i = 0; i < 3; ++i)
        {
            relative[i] /= runs;
            absolute[i] /= runs;
        }
        
        printf("| %4d | %e | %e | %e | absolute error  |\n", bw, absolute[0], absolute[1], absolute[2]);
        printf("|      | %e | %e | %e | relative error  |\n", relative[0], relative[1], relative[2]);
        printf("+------+--------------+--------------+--------------+-----------------+\n");
    }
    
    return 0;
}