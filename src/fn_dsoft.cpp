//
//  fn_dsoft.cpp
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

#include <pfsoft>

PFSOFT_NAMESPACE(FourierTransforms)

/*!
 * @brief           The DSOFT (<b>S0</b>(3) <b>F</b>ourier <b>T</b>ransform)
 *                  describes the FFT on the rotation group \f$\mathcal{SO}(3)\f$
 * @details         The method to compute the DSOFT on the rotation group is described in
 *                  detail in the paper 'FFTs on the Rotation Group' written by Peter J.
 *                  Kostelec and Daniel N. Rockmore. The implementation that is underneath
 *                  this function is
 *                  \f[
 *                      \hat{f}^l_{M,M'} = \frac{\pi}{(2B)^2}\sum\limits_{k = 0}^{2B-1}w_B(k)
 *                          \tilde{d}^l_{M,M'}(\beta_k)\sum\limits_{j2 = 0}^{2B-1}e^{iM'\gamma_{j_2}}
 *                          \sum\limits_{j_1 = 0}^{2B-1}e^{iM\alpha_{j_1}}f(\alpha_{j_1},\beta_k,\gamma_{j_2})
 *                  \f]
 *                  where \f$B\f$ is the bandwidth of function \f$f(\alpha_{j_1},\beta_k,\gamma_{j_2})\f$.
 *                  The number of cofficients \f$\hat{f}^l_{M,M'}\f$ can be calculated by
 *                  \f[
 *                      |\{\hat{f}^l_{M,M'}\}_{l\in\{0,\dots,B-1\}, M,M'\in\{-l,\dots,l\}}| = \sum\limits_{i = 0}^{B-1}(2 * i + 1)^2
 *                  \f]
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
 * @param[in]       sample A discrete sample of function \f$f\f$ which has the
 *                  dimension of \f$2B\times 2B\times 2B\f$.
 * @param[out]      fc A Fourier coefficent managment container with capacaty for
 *                  all Fourier coefficients of \f$f\f$.
 *
 * @sa              DWT::quadrature_weights
 * @sa              DWT::wigner_d_matrix
 * @sa              DSOFTFourierCoefficients
 * @sa              complex
 * @sa              matrix
 * @sa              grid3D
 *
 * @ingroup         FourierTransforms
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            14.05.2015
 *
 * @since           0.0.1
 */
void DSOFT(grid3D< complex< double > > sample, DSOFTFourierCoefficients& fc, int threads)
{
    /*****************************************************************
     ** Check parameters                                            **
     *****************************************************************/
    // Check if the grid has same size in each dimension
    pfsoft_cond_w_ret(sample.rows != sample.cols || sample.rows != sample.lays, "%s", "all DSOFT sample grid dimensions should be equal.");
    
    // Check if grid has odd dimensions
    pfsoft_cond_w_ret(sample.rows & 1, "%s", "DSOFT sample grid dimensions are not even.");
    
    // Extract bandwidth
    const int bandwidth = static_cast< int >(sample.cols / 2);
    
    // precompute the double bandwidth
    const int bw2 = 2 * bandwidth;
    
    // Check if Fourier coefficients container dimension matches sample dimension
    pfsoft_cond_w_ret(bandwidth != fc.bandwidth, "%s", "DSOFT Fourier coefficients container bandwidth does not match to sample grid bandwidth.");
    
    // print warinings for serial implementation
    #ifndef _OPENMP
    pfsoft_cond_w(threads != 1, "%s", "compiler does not support OpenMP. Changing the number of threads for the DSOFT has no effect.");
    #endif
    
    /*****************************************************************
     ** FFT2 transform layers of sample grid for fixed k            **
     *****************************************************************/
    sample.layer_wise_DFT2();
    
    /*****************************************************************
     ** M = 0, M' = 0                                               **
     *****************************************************************/
    vector< long double > weights(2 * bandwidth);
    DWT::quadrature_weights< long double >(weights);
    
    matrix< long double > dw(bandwidth, 2 * bandwidth);
    DWT::weighted_wigner_d_matrix(dw, bandwidth, 0, 0, weights);
    dw *= -1;
    
    vector< complex< double > > s(bw2, vector< complex< double > >::COLUMN);
    
    // defining norm factor
    complex< double > norm(constants< double >::pi / (bandwidth * bw2), 0);
    
    // defining needed indices
    int MMp, M, Mp;
    vector< complex< double > >::iterator e;
    matrix< long double >::lin_iterator m;
    
    // DWT for M = 0, M' = 0
    for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(0, 0, e - s.begin());                            }
    vector< complex< double > > sh = dw * s;
    for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), 0, 0) = norm * *e;            }
    
    /*****************************************************************
     ** Iterate over all combinations of M and M'                   **
     *****************************************************************/
    #pragma omp parallel default(shared) if(bandwidth >= DSOFT_THRESHOLD) num_threads(threads)
    {
        
        #pragma omp for private(M, e) firstprivate(dw, s, sh) schedule(dynamic) nowait
        for (M = 1; M < bandwidth; ++M)
        {
            access::rw(dw.rows) = bandwidth - M;
            DWT::weighted_wigner_d_matrix(dw, bandwidth, M, 0, weights);
            dw *= -1;
            
            /*****************************************************************
             ** Make use of symmetries                                      **
             *****************************************************************/
            // case f_{M,0}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(0, M, e - s.begin());                    }
            sh = dw * s;
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), M, 0) = norm * *e;    }
            
            // case f_{0,M}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(M, 0, e - s.begin());                    }
            sh = dw * s;
            if  (M & 1)                              { sh *= -1;                                            }
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), 0, M) = norm * *e;    }
            
            // case f_{-M,0}
            fliplr(dw);
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(0, bw2 - M, e - s.begin());              }
            sh = dw * s;
            if (M & 1)  // if M is odd
            {
                for (e = sh.begin(); e < sh.end(); e += 2)     { *e *= -1;                                  }
            }
            else        // if M is even
            {
                for (e = sh.begin() + 1; e < sh.end(); e += 2) { *e *= -1;                                  }
            }
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), -M, 0) = norm * *e;   }
            
            // case f_{0,-M}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(bw2 - M, 0, e - s.begin());              }
            sh = dw * s;
            for (e = sh.begin() + 1; e < sh.end(); e += 2) { *e *= -1;                                      }
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), 0, -M) = norm * *e;   }
            
            // get new wigner matrix
            DWT::weighted_wigner_d_matrix(dw, bandwidth, M, M, weights);
            dw *= -1;
            
            // case f_{M, M}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(M, M, e - s.begin());                    }
            sh = dw * s;
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), M, M) = norm * *e;    }
            
            // case f_{-M, -M}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(bw2 - M, bw2 - M, e - s.begin());        }
            sh = dw * s;
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), -M, -M) = norm * *e;  }
            
            // Modify dw for the last two cases. flip matrix from left to right and negate signs of
            // every second row with odd row indices.
            fliplr_ne2ndorow(dw);
            
            // A little arithmetic error is occuring in the following calculation... I do not exactly know why
            // case f_{M, -M}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(bw2 - M, M, e - s.begin());              }
            sh = dw * s;
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), M, -M) = norm * *e;   }
            
            // case f_{-M, M}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(M, bw2 - M, e - s.begin());              }
            sh = dw * s;
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), -M, M) = norm * *e;   }
        }
        
        // Fused two loops per hand
        //
        // for (M = 1; M < bandwidth; ++M)
        //     for (Mp = 1; Mp < M; ++Mp)
        //
        // which now is equivalent to the following loop
        #pragma omp for private(MMp, M, Mp, e) firstprivate(dw, s, sh) schedule(dynamic) nowait
        for (MMp = 0; MMp < (bandwidth - 2) * (bandwidth - 1) / 2; ++MMp)
        {
            // reconstructing nested loop indices
            int i = MMp / (bandwidth - 1) + 1;
            int j = MMp % (bandwidth - 1) + 1;
            
            // get M and M'
            M  = j > i ? bandwidth - i : i + 1;
            Mp = j > i ? bandwidth - j : j    ;
            
            // get new wigner d-matrix
            access::rw(dw.rows) = bandwidth - std::max(abs(M), abs(Mp));
            DWT::weighted_wigner_d_matrix(dw, bandwidth, M, Mp, weights);
            
            // case f_{M, Mp}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(Mp, M, e - s.begin());                   }
            sh  = dw * s;
            sh *= -1;
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), M, Mp) = norm * *e;   }
            
            // case f_{Mp, M}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(M, Mp, e - s.begin());                   }
            sh = dw * s;
            if  (!((M - Mp) & 1))                    { sh *= -1;                                            }
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), Mp, M) = norm * *e;   }
            
            // case f_{-M, -Mp}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(bw2 - Mp, bw2 - M, e - s.begin());       }
            sh = dw * s;
            if  (!((M - Mp) & 1))                    { sh *= -1;                                            }
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), -M, -Mp) = norm * *e; }
            
            // case f_{-Mp, -M}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(bw2 - M, bw2 - Mp, e - s.begin());       }
            sh  = dw * s;
            sh *= -1;
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), -Mp, -M) = norm * *e; }
            
            // modify wigner d-matrix for next four cases. This just works because the weight
            // function is also symmetric like the wigner-d matrix. flip left-right the dw
            // matrix and negate each even value with even row index.
            fliplr_ne2nderow(dw);
            
            // case f_{Mp, -M}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(bw2 - M, Mp, e - s.begin());             }
            sh = dw * s;
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), Mp, -M) = norm * *e;  }
            
            // case f_{M, -Mp}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(bw2 - Mp, M, e - s.begin());             }
            sh = dw * s;
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), M, -Mp) = norm * *e;  }
            
            // alter signs
            if ((M - Mp) & 1)
            {
                for (m = dw.begin(); m != dw.begin() + dw.rows * dw.cols; ++m) { *m *= -1;                  }
            }
            
            // case f_{-Mp, M}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(M, bw2 - Mp, e - s.begin());             }
            sh = dw * s;
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), -Mp, M) = norm * *e;  }
            
            // case f_{-M, Mp}
            for (e = s.begin() ; e != s.end() ; ++e) { *e = sample(Mp, bw2 - M, e - s.begin());             }
            sh = dw * s;
            for (e = sh.begin(); e != sh.end(); ++e) { fc(bandwidth - (sh.end() - e), -M, Mp) = norm * *e;  }
        }
    }
}

PFSOFT_NAMESPACE_END
