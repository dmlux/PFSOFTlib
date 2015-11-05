//
//  fn_idsoft.cpp
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

PDSOFT_NAMESPACE(FourierTransforms)

/*!
 * @brief           The inverse DSOFT (<b>S0</b>(3) <b>F</b>ourier <b>T</b>ransform)
 *                  describes the inverse FFT on the rotation group \f$\mathcal{SO}(3)\f$
 * @details         The method to compute the inverse DSOFT on the rotation group is described
 *                  in detail in the paper 'FFTs on the Rotation Group' written by Peter J.
 *                  Kostelec and Daniel N. Rockmore. The implementation that is underneath
 *                  this function is
 *                  \f[
 *                      f(\alpha,\beta,\gamma) = \sum\limits_{J\geq 0}\sum\limits^J_{M=-J}\sum\limits^J_{M'=-J}
 *                          \hat{f}^J_{MM'}\tilde{D}^J_{MM'}(\alpha,\beta,\gamma)
 *                  \f]
 *                  where \f$\hat{f}^J_{MM'}\f$ is the Fourier coefficient of degree
 *                  \f$J\f$ and orders \f$M,\;M'\f$. For more detailed information about the
 *                  Fourier coefficients read the documentation of FourierTransforms::DSOFT.
 *                  The implementation itself uses the symmetry properties of the wigner
 *                  d-function to reduce the number of wigner d-function evaluations. The
 *                  symmetries that are used are
 *                  \f{eqnarray*}{
 *                      d^{J}_{MM'}(\beta) &=& (-1)^{M-M'}d^J_{-M-M'}(\beta)\\
 *                      &=& (-1)^{M-M'}d^J_{M'M}(\beta)\\
 *                      &=& d^J_{-M'-M}(\beta)\\
 *                      &=& (-1)^{J-M'}d^J_{-MM'}(\pi-\beta)\\
 *                      &=& (-1)^{J+M}d^J_{M-M'}(\pi-\beta)\\
 *                      &=& (-1)^{J-M'}d^J_{-M'M}(\pi-\beta)\\
 *                      &=& (-1)^{J+M}d^J_{M'-M}(\pi-\beta)
 *                  \f}
 *
 * @param[in]       fc A Fourier coefficent managment container with all Fourier coefficients
 *                  of the DSOFT.
 * @param[out]      synthesis The synthesized sample for the given Fourier coefficients.
 *
 * @sa              DWT::wigner_d_matrix
 * @sa              DSOFTFourierCoefficients
 * @sa              FourierTransforms::DSOFT
 * @sa              complex
 * @sa              matrix
 * @sa              grid3D
 *
 * @since           0.0.1
 *
 * @ingroup         FourierTransforms
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            23.05.2015
 */
void IDSOFT(const DSOFTFourierCoefficients& fc, grid3D< complex< double > >& synthesis, int threads)
{
    /*****************************************************************
     ** Check parameters                                            **
     *****************************************************************/
    // Check if the grid has same size in each dimension
    if (synthesis.rows != synthesis.cols || synthesis.rows != synthesis.lays)
    {
        pdsoft_warning("%s", "all ISOFT synthesis grid dimensions should be equal.");
        return;
    }
    
    // Check if grid has odd dimensions
    if (synthesis.rows & 1)
    {
        pdsoft_warning("%s", "ISOFT synthesis grid dimensions are not even.");
        return;
    }
    
    // Extract bandwidth
    int bandwidth = static_cast< int >(synthesis.cols / 2);
    
    // precompute the double bandwidth
    const int bw2 = 2 * bandwidth;
    
    // Check if Fourier coefficients container dimension matches sample dimension
    if (bandwidth != fc.bandwidth)
    {
        pdsoft_warning("%s", "ISOFT Fourier coefficients container bandwidth does not match to synthesis grid bandwidth.");
        return;
    }
    
    // print warinings for serial implementation
    #ifndef _OPENMP
    if (threads != 1)
    {
        pdsoft_warning("%s", "compiler does not support OpenMP. Changing the number of threads for the ISOFT has no effect.");
    }
    #endif
    
    /*****************************************************************
     ** M = 0, M' = 0                                               **
     *****************************************************************/
    matrix< double > d(bandwidth, 2 * bandwidth);
    DWT::wigner_d_matrix< double >(d, bandwidth, 0, 0);
    
    d *= -1;
    d.transpose();
    
    vector< complex< double > > sh(d.cols, vector< complex< double > >::COLUMN);
    
    // defining norm factor
    complex< double > norm((bandwidth * bw2) / constants< double >::pi, 0);
    
    // defining needed indices
    int MMp, M, Mp;
    
    // defining type for following iterations
    typedef const complex< double >* cx_it;
    
    // inverse DWT for M = 0, M' = 0
    for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), 0, 0);  }
    vector< complex< double > > s = d * sh;
    for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(0, 0, e - s.mem) = *e;                                       }
    
    /*****************************************************************
     ** Iterate over all combinations of M and M'                   **
     *****************************************************************/
    #pragma omp parallel default(shared) if(bandwidth >= DSOFT_THRESHOLD) num_threads(threads)
    {
        
        #pragma omp for private(M, d, s, sh) schedule(dynamic) nowait
        for (M = 1; M < bandwidth; ++M)
        {
            d  = matrix< double >(bandwidth - M, 2 * bandwidth);
            DWT::wigner_d_matrix< double >(d, bandwidth, M, 0);
            
            d *= -1;
            d.transpose();
            
            sh = vector< complex< double > >(d.cols, vector< complex< double > >::COLUMN);
            
            /*****************************************************************
             ** Make use of symmetries                                      **
             *****************************************************************/
            // case f_{M,0}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), M, 0);      }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(0, M, e - s.mem) = *e;                                           }
            
            // case f_{0,M}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), 0, M);      }
            if  (M & 1) { s = d * (sh * -1); } else                { s = d * sh;                                                                }
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(M, 0, e - s.mem) = *e;                                           }
            
            flipud(d);
            
            // case f_{-M,0}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), -M, 0);     }
            if (M & 1)
            {
                for (cx_it e = sh.mem; e < sh.mem + sh.size; e+=2)     { access::rw(*e) *= -1;                                                  }
            }
            else
            {
                for (cx_it e = sh.mem + 1; e < sh.mem + sh.size; e+=2) { access::rw(*e) *= -1;                                                  }
            }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(0, bw2 - M, e - s.mem) = *e;                                     }
            
            // case f_{0,-M}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), 0, -M);     }
            for (cx_it e = sh.mem + 1; e < sh.mem + sh.size; e+=2) { access::rw(*e) *= -1;                                                      }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(bw2 - M, 0, e - s.mem) = *e;                                     }
            
            // get new wigner matrix
            d  = matrix< double >(bandwidth - M, 2 * bandwidth);
            DWT::wigner_d_matrix< double >(d, bandwidth, M, M);
            
            d *= -1;
            d.transpose();
            
            // case f_{M,M}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), M, M);      }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(M, M, e - s.mem) = *e;                                           }
            
            // case f_{-M,-M}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), -M, -M);    }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(bw2 - M, bw2 - M, e - s.mem) = *e;                               }
            
            // Modify dw for the last two cases. flip matrix from left to right and negate every
            // second row with odd row indices.
            flipud_ne2ndocol(d);
            
            // An little arithmetic error is occuring in the following calculation... I do not exactly know why...
            // case f_{M,-M}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), M, -M);     }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(bw2 - M, M, e - s.mem) = *e;                                     }
            
            // case f_{-M,M}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), -M, M);     }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(M, bw2 - M, e - s.mem) = *e;                                     }
        }
        
        // Fused two loops per hand
        //
        // for (M = 1; M < bandwidth; ++M)
        //     for (Mp = 1; Mp < M; ++Mp)
        //
        // which now is equivalent to the following loop
        #pragma omp for private(MMp, M, Mp, d, s, sh) schedule(dynamic) nowait
        for (MMp = 0; MMp < (bandwidth - 2) * (bandwidth - 1) / 2; ++MMp)
        {
            // reconstructing indices of the two nested for loops
            int i = MMp / (bandwidth - 1) + 1;
            int j = MMp % (bandwidth - 1) + 1;
            
            // get M and M'
            M  = j > i ? bandwidth - i : i + 1;
            Mp = j > i ? bandwidth - j : j    ;
            
            // get new wigner d-matrix
            d  = matrix< double >(bandwidth - std::max(abs(M), abs(Mp)), 2 * bandwidth);
            DWT::wigner_d_matrix< double >(d, bandwidth, M, Mp);
            d.transpose();
            
            sh = vector< complex< double > >(d.cols, vector< complex< double > >::COLUMN);
            
            // case f_{M,Mp}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), M, Mp);     }
            for (cx_it e = sh.mem; e < sh.mem + sh.size; ++e)      { access::rw(*e) *= -1;                                                      }
            s  = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(Mp, M, e - s.mem) = *e;                                          }
            
            // case f_{Mp,M}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), Mp, M);     }
            if  (!((M - Mp) & 1))                                  { sh *= -1;                                                                  }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(M, Mp, e - s.mem) = *e;                                          }
            
            // case f_{-M,-Mp}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), -M, -Mp);   }
            if  (!((M - Mp) & 1))                                  { sh *= -1;                                                                  }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(bw2 - Mp, bw2 - M, e - s.mem) = *e;                              }
            
            // case f_{-Mp,-M}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), -Mp, -M);   }
            for (cx_it e = sh.mem; e < sh.mem + sh.size; ++e)      { access::rw(*e) *= -1;                                                      }
            s  = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(bw2 - M, bw2 - Mp, e - s.mem) = *e;                              }
            
            // modify wigner d-matrix for next four cases. This just works because the weight
            // function is also symmetric like the wigner-d matrix. flip up-dow the d
            // matrix and negate every second column with even row index.
            flipud_ne2ndecol(d);
            
            // case f_{Mp,-M}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), Mp, -M);    }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(bw2 - M, Mp, e - s.mem) = *e;                                    }
            
            // case f_{M,-Mp}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), M, -Mp);    }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(bw2 - Mp, M, e - s.mem) = *e;                                    }
            
            // alter signs
            if ((M - Mp) & 1)
            {
                for (const double* e = d.mem; e < d.mem + d.rows * d.cols; ++e) { access::rw(*e) *= -1;                                         }
            }
            
            // case f_{-Mp,M}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), -Mp, M);    }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(M, bw2 - Mp, e - s.mem) = *e;                                    }
            
            // case f_{-M,Mp}
            for (cx_it e = sh.mem + sh.size - 1; e >= sh.mem; --e) { access::rw(*e) = norm * fc(bandwidth - (sh.mem + sh.size - e), -M, Mp);    }
            s = d * sh;
            for (cx_it e = s.mem + bw2 - 1; e >= s.mem; --e)       { synthesis(Mp, bw2 - M, e - s.mem) = *e;                                    }
        }
    }
    
    /*****************************************************************
     ** IFFT2 transform layers of input sample grid for fixed k     **
     *****************************************************************/
    synthesis.layer_wise_IDFT2(complex< double > (1. / (4. * bandwidth * bandwidth), 0));
}

PDSOFT_NAMESPACE_END
