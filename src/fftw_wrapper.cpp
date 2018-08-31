//
//  fftw_wrapper.cpp
//  PFSOFTlib
//
//   Created by Denis-Michael Lux on 05. November 2015.
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

/*- FFTW3 -*/
/*- for Fourier Analysis purposes -*/
#include <fftw3.h>

/*- including function wrapper    -*/
#include <pfsoft>

/*!
 * @brief           The namespace containg all includes of datatypes, functions and 
 *                  other namespaces for mathmatical computation from the **UZLMath**
 *                  library
 * @details         All datatypes and functions that are created for the library
 *                  'UZLMath' are encapsulated in this namespace.
 *
 * @since           0.0.1
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            14.05.2015
 */
PFSOFT_BEGIN

/*- Wrapper for needed library functions -*/
extern "C"
{
    void uzl_fftw_layer_wise_DFT2_grid3D(int cols, int rows, int lays, double* arr, int threads)
    {
        #ifdef _OPENMP
        
        // define indices
        int i;
        
        // storage plans
        fftw_plan plans[lays];
        
        // create plans
        // Creating plans and destroying them is not really thread safe
        // therefore the plans should be created and destoryed serially
        for (i = 0; i < lays; ++i)
        {
            // get correct layer
            fftw_complex* layer = (fftw_complex*)arr + i * rows * cols;
            
            // create plan to execute an FFT
            plans[i] = fftw_plan_dft_2d(cols, rows, layer, layer, FFTW_FORWARD, FFTW_ESTIMATE);
        }
        
        // execute plans in parallel
        #pragma omp parallel for private(i) shared(lays, plans) schedule(dynamic) num_threads(threads)
        for (i = 0; i < lays; ++i)
        {
            // execute FFT2 plan
            fftw_execute(plans[i]);
        }
        
        for (i = 0; i < lays; ++i)
        {
            fftw_destroy_plan(plans[i]);
        }
        
        // remove all additional allocated objects needed for the FFT
        fftw_cleanup();
        
        #else
        
        // define parameters for many dft
        int rank    = 2;            // Dimension of FFT's
        int n[]     = {cols, rows}; // Dimension of matrix for each FFT
        int howmany = lays;         // How many FFT's should be executed?
        int idist   = n[0] * n[1];  // Length of memory space for each matrix
        int odist   = idist;        // ... for output same, but is in-place...
        int istride = 1;            // 1 because matrices are contigous in memory
        int ostride = 1;            // ... for output same, but is still in-place...
        int* inembed= n;
        int* onembed= n;
        
        // create plan to execute an layerwise FFT2
        fftw_plan lay_wise_fft2 = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)arr, inembed, istride, idist, (fftw_complex*)arr, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
        
        // execute layer-wise FFT2
        fftw_execute(lay_wise_fft2);
        
        // free allocated memory
        fftw_destroy_plan(lay_wise_fft2);
        fftw_cleanup();
        
        #endif
    }
    
    void uzl_fftw_layer_wise_IDFT2_grid3D(int cols, int rows, int lays, double* arr, int threads)
    {
        #ifdef _OPENMP
        
        // define indices
        int i;
        
        // storage plans
        fftw_plan plans[lays];
        
        // create plans
        // Creating plans and destroying them is not really thread safe
        // therefore the plans should be created and destoryed serially
        for (i = 0; i < lays; ++i)
        {
            // get correct layer
            fftw_complex* layer = (fftw_complex*)arr + i * rows * cols;
            
            // create plan to execute an FFT
            plans[i] = fftw_plan_dft_2d(cols, rows, layer, layer, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
        
        // execute plans in parallel
        #pragma omp parallel for private(i) shared(lays, plans) schedule(dynamic) num_threads(threads)
        for (i = 0; i < lays; ++i)
        {
            // execute FFT2 plan
            fftw_execute(plans[i]);
        }
        
        for (i = 0; i < lays; ++i)
        {
            fftw_destroy_plan(plans[i]);
        }
        
        // remove all additional allocated objects needed for the FFT
        fftw_cleanup();
        
        #else
        
        // define parameters for many dft
        int rank    = 2;            // Dimension for IFFT's
        int n[]     = {cols, rows}; // Dimension of matrix for each IFFT
        int howmany = lays;         // How many IFFT's should be executed?
        int idist   = n[0] * n[1];  // Length of memory space for each matrix
        int odist   = idist;        // ... for output same, but is in-place...
        int istride = 1;            // 1 because matrices are contigous in memory
        int ostride = 1;            // ... for output same, but is still in-place...
        int* inembed= n;
        int* onembed= n;
        
        // create plan to execute an layerwise IFFT2
        fftw_plan lay_wise_ifft2 = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)arr, inembed, istride, idist, (fftw_complex*)arr, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
        
        // execute layer-wise IFFT2
        fftw_execute(lay_wise_ifft2);
        
        // free allocated memory
        fftw_destroy_plan(lay_wise_ifft2);
        fftw_cleanup();
        
        #endif
    }
}

PFSOFT_END
