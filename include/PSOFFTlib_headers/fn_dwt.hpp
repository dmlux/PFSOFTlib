//
//  fn_dwt.hpp
//  PSOFFTlib
//
//   Created by Denis-Michael Lux on 05. November 2015.
//
//   This file is part of PSOFFTlib.
//
//   PSOFFTlib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   PSOFFTlib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with PSOFFTlib.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef PSOFFTlib_fn_dwt_hpp
#define PSOFFTlib_fn_dwt_hpp

PSOFFT_NAMESPACE(DWT)
    
/*- For more information/implementation details see fn_dwt.cpp file! -*/

/*!
 * @brief       Computes the **quadrature weights** as a vector that is used for the DWT
 *              which is necessary for a bandwidth \f$B\f$ transform.
 * @details     The weights vector is a vector that contains the needed weights for the
 *              data vector. The elements can be expressed by
 *              \f{eqnarray*}{
 *                  w_B(j) = \frac{2}{B}\sin\left(\frac{\pi(2j+1)}{4B}\right)
 *                      \sum\limits_{k = 0}^{B-1}\frac{1}{2k + 1}
 *                      \sin\left((2j+1)(2k+1)\frac{\pi}{4B}\right)
 *              \f}
 *              where \f$B\f$ is the bandlimit that is given and \f$0\leq j\leq 2B-1\f$.
 *              The dimension of this matrix is \f$2B\times 2B\f$
 *
 * @param[in]   bandwidth The given bandwidth
 * @return      A vector containing the quadrature weights that can be used to compute the DWT
 *              \f{eqnarray*}{
 *                  \begingroup
 *                  \renewcommand*{\arraystretch}{1.5}
 *                  W = \left(\begin{array}{c}
 *                      w_B(0)\\
 *                      w_B(1)\\
 *                      \vdots\\
 *                      w_B(2B-1)
 *                  \end{array}\right)
 *                  \endgroup
 *              \f}
 *
 * @since       0.0.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        03.05.15
 */
template< typename T >
inline
void_number_type< T > quadrature_weights(vector< T >& vec)
{
    typedef T pod_type;
    psofft_cond_w_ret(vec.size & 1, "%s", "uneven vector length in DWT::quadrature_weights. ");
    
    int i, k, bandwidth = vec.size / 2;
    for (i = 0; i < bandwidth; ++i)
    {
        pod_type wi  = 2.0 / bandwidth * sin(constants< T >::pi * (2.0 * i + 1.0)/(4.0 * bandwidth));
        pod_type sum = 0;
        for (k = 0; k < bandwidth; ++k)
        {
            sum += 1.0 / (2.0 * k + 1.0) * sin((2.0 * i + 1.0) * (2.0 * k + 1.0) * constants< T >::pi / (4.0 * bandwidth));
        }
        
        wi                        *= sum;
        vec[i]                     = wi;
        vec[2 * bandwidth - 1 - i] = wi;
    }
}

/*!
 * @brief       The Wigner d-matrix where the weights are calculated onto the matrix values.
 * @details     Calculates \f$d\cdot w\f$ where \f$d\f$ is matrix containing wigner
 *              d-Function values on each entry and \f$w\f$ is diagonal matrix containing
 *              the quadrature weights on the diagonal.
 *
 * @param[in]   bandwidth The given bandwidth.
 * @param[in]   M The order \f$M\f$ of \f$d^J_{MM'}\f$.
 * @param[in]   Mp The order \f$M'\f$ of \f$d^J_{MM'}\f$.
 * @param[in]   weights A vector containing the quadrature weights.
 * @return      A matix containing weighted values of Wigner d-function.
 *
 * @sa          wigner::wiger_d_matrix
 * @sa          wigner::weight_matrix
 *
 * @since       0.0.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        03.05.15
 */
template< typename T >
inline
void_number_type< T > weighted_wigner_d_matrix(matrix< T >& wig, const int& bandwidth, const int& M, const int& Mp, const vector< T >& weights)
{
    typedef T pod_type;
    
    // Definition of used indices and the matrix that will be returned
    int i, j, minJ = std::max(abs(M), abs(Mp));
    
    psofft_cond_w_ret(wig.rows != bandwidth - minJ || wig.cols != 2 * bandwidth, "%s", "dimension mismatch between input matrix and function arguments in DWT::weighted_wigner_d_matrix.");
    
    // Compute root coefficient for the base case
    pod_type normFactor  = sqrt((2.0 * minJ + 1.0)/2.0);
    for (i = 0 ; i < minJ - std::min(abs(M), abs(Mp)) ; ++i)
    {
        normFactor *= sqrt((2.0 * minJ - i) / (i + 1.0));
    }
    
    // Sin sign for the recurrence base case
    pod_type sinSign = (minJ == abs(M) && M >= 0 && (minJ - Mp) & 1 ? 1       : -1);
    sinSign          = (minJ != abs(M) && Mp < 0 && (minJ - Mp) & 1 ? sinSign : -1);
    
    // Powers
    pod_type cosPower, sinPower;
    if (minJ == abs(M) && M >= 0)
    {
        cosPower = minJ + Mp;
        sinPower = minJ - Mp;
    }
    else if (minJ == abs(M))
    {
        cosPower = minJ - Mp;
        sinPower = minJ + Mp;
    }
    else if (Mp >= 0)
    {
        cosPower = minJ + M;
        sinPower = minJ - M;
    }
    else
    {
        cosPower = minJ - M;
        sinPower = minJ + M;
    }
    
    // Base cases and filling matrix with values
    pod_type cosBeta[2 * bandwidth];
    for (i = 0 ; i < 2 * bandwidth; ++i)
    {
        // Getting sin and cos values for the power operator
        pod_type sinHalfBeta = sin(0.5 * ((2.0 * i + 1.0) * constants< T >::pi) / (4.0 * bandwidth));
        pod_type cosHalfBeta = cos(0.5 * ((2.0 * i + 1.0) * constants< T >::pi) / (4.0 * bandwidth));
        
        // Store cosine values for reuse in recurrence loop
        cosBeta[i] = cos(((2.0 * i + 1.0) * constants< T >::pi) / (4.0 * bandwidth));
        
        // Computing base wigners. filling the first row in matrix with those values
        wig(0, i)  = normFactor * sinSign * pow(sinHalfBeta, sinPower) * pow(cosHalfBeta, cosPower) * weights[i];
    }
    
    // Filling wigner matrix with values. Starting with second row and
    // iterate to last row with index B - 1
    for(i = 0 ; i < bandwidth - minJ - 1; ++i)
    {
        // Recurrence coefficients
        pod_type c1   = 0;
        
        // Index for wigner function
        pod_type idx  = minJ + i;
        
        // Terms in recurrence
        pod_type norm = sqrt((2.0 * idx + 3.0) / (2.0 * idx + 1.0));
        pod_type nom  = (idx + 1.0) * (2. * idx + 1.0);
        pod_type den  = 1.0 / sqrt(((idx + 1) * (idx + 1) - M*M) * ((idx + 1) * (idx + 1) - Mp*Mp));
        
        // Fractions
        pod_type f1   = norm * nom * den;
        pod_type f2   = 0;
        
        // Correcting undefined values from division by zero
        if (minJ + i != 0)
        {
            pod_type t1 = sqrt((2.0 * idx + 3.0)/(2.0 * idx - 1.0) ) * (idx + 1.0)/idx ;
            pod_type t2 = sqrt((idx*idx - M*M) * (idx*idx - Mp*Mp));
            
            c1   = -t1 * t2 * den;
            f2   = -M*Mp / (idx * (idx + 1.));
        }
        
        //  Filling matrix with next recurrence step value
        for (j = 0; j < 2*bandwidth; ++j)
        {
            wig(i + 1, j) = c1 * (i == 0 ? 0 : wig(i - 1, j)) + wig(i, j) * f1 * (f2 + cosBeta[j]);
        }
    }
}

/*!
 * @brief       Generates the Wigner d-matrix which is a matrix that
 *              contains the results of the \f$L^2\f$-normalized Wigner
 *              d-function on each entry.
 * @details     The dimension of this matrix is \f$(B-J+1)\times 2B\f$ where
 *              \f$B\f$ denotes a given bandwidth.
 *
 * @param[in]   bandwidth The given bandwidth
 * @param[in]   M The order \f$M\f$ for the \f$L^2\f$-normalized Wigner d-function
 * @param[in]   Mp The order \f$M'\f$ for the \f$L^2\f$-normalized Wigner d-function
 * @return      The resulting matrix looks as follows
 *              \f{eqnarray*}{
 *                  \begingroup
 *                  \renewcommand*{\arraystretch}{1.5}
 *                  D = \left(\begin{array}{c c c c}
 *                      \tilde{d}^J_{M,M'}(\beta_0)      & \tilde{d}^J_{M,M'}(\beta_1)
 *                      & \cdots                        & \tilde{d}^J_{M,M'}(\beta_{2B-1})\\
 *                      \tilde{d}^{J+1}_{M,M'}(\beta_0)  & \tilde{d}^{J+1}_{M,M'}(\beta_1)
 *                      & \cdots                        & \tilde{d}^{J+1}_{M,M'}(\beta_{2B-1})\\
 *                      \vdots & \vdots & \ddots & \vdots\\
 *                      \tilde{d}^{B-1}_{M,M'}(\beta_0)  & \tilde{d}^{B-1}_{M,M'}(\beta_1)
 *                      & \cdots                        & \tilde{d}^{B-1}_{M,M'}(\beta_{2B-1})
 *                  \end{array}\right)
 *                  \endgroup
 *              \f}
 *              with \f$\beta_k = \frac{\pi(2k + 1)}{4B}\f$
 *
 * @since       0.0.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        03.05.15
 *
 * @see         wigner::wigner_d
 * @see         wigner::wigner_d_l2normalized
 */
template< typename T >
inline
void_number_type< T > wigner_d_matrix(matrix< T >& wig, const int& bandwidth, const int& M, const int& Mp)
{
    typedef T pod_type;
    
    // Definition of used indices and the matrix that will be returned
    int i, j, minJ = std::max(abs(M), abs(Mp));
    
    psofft_cond_w_ret(wig.rows != bandwidth - minJ || wig.cols != 2 * bandwidth, "%s", "dimension mismatch between input matrix and function arguments in DWT::weighted_wigner_d_matrix.");
    
    // Compute root coefficient for the base case
    pod_type normFactor  = sqrt((2.0 * minJ + 1.0)/2.0);
    for (i = 0 ; i < minJ - std::min(abs(M), abs(Mp)) ; ++i)
    {
        normFactor *= sqrt((2.0 * minJ - i) / (i + 1.0));
    }
    
    // Sin sign for the recurrence base case
    pod_type sinSign = (minJ == abs(M) && M >= 0 && (minJ - Mp) & 1 ? 1 : -1      );
    sinSign          = (minJ != abs(M) && Mp < 0 && (minJ - Mp) & 1 ? sinSign : -1);
    
    // Powers
    pod_type cosPower, sinPower;
    if (minJ == abs(M) && M >= 0)
    {
        cosPower = minJ + Mp;
        sinPower = minJ - Mp;
    }
    else if (minJ == abs(M))
    {
        cosPower = minJ - Mp;
        sinPower = minJ + Mp;
    }
    else if (Mp >= 0)
    {
        cosPower = minJ + M;
        sinPower = minJ - M;
    }
    else
    {
        cosPower = minJ - M;
        sinPower = minJ + M;
    }
    
    // Base cases and filling matrix with values
    pod_type cosBeta[2 * bandwidth];
    for (i = 0 ; i < 2 * bandwidth; ++i)
    {
        // Getting sin and cos values for the power operator
        pod_type sinHalfBeta = sin(0.5 * ((2.0 * i + 1.0) * constants< T >::pi) / (4.0 * bandwidth));
        pod_type cosHalfBeta = cos(0.5 * ((2.0 * i + 1.0) * constants< T >::pi) / (4.0 * bandwidth));
        
        // Store cosine values for reuse in recurrence loop
        cosBeta[i] = cos(((2.0 * i + 1.0) * constants< T >::pi) / (4.0 * bandwidth));
        
        // Computing base wigners
        wig(0, i)  = normFactor * sinSign * pow(sinHalfBeta, sinPower) * pow(cosHalfBeta, cosPower);
    }
    
    // Filling wigner matrix with values. Starting with second row and
    // iterate to last row with index B - 1
    for(i = 0 ; i < bandwidth - minJ - 1; ++i)
    {
        // Recurrence coefficients
        pod_type c1   = 0;
        
        // Index for wigner function
        pod_type idx  = minJ + i;
        
        // Terms in recurrence
        pod_type norm = sqrt((2.0 * idx + 3.0) / (2.0 * idx + 1.0));
        pod_type nom  = (idx + 1.0) * (2.0 * idx + 1.0);
        pod_type den  = 1.0 / sqrt(((idx + 1) * (idx + 1) - M*M) * ((idx + 1) * (idx + 1) - Mp*Mp));
        
        // Fractions
        pod_type f1   = norm * nom * den;
        pod_type f2   = 0;
        
        // Correcting undefined values from division by zero
        if (minJ + i != 0)
        {
            T t1 = sqrt((2.0 * idx + 3.0)/(2.0 * idx - 1.0) ) * (idx + 1.0)/idx ;
            T t2 = sqrt((idx*idx - M*M) * (idx*idx - Mp*Mp));
            
            c1   = -t1 * t2 * den;
            f2   = -M*Mp / (idx * (idx + 1.0));
        }
        
        //  Filling matrix with next recurrence step value
        for (j = 0; j < 2 * bandwidth; ++j)
        {
            wig(i + 1, j) = c1 * (i == 0 ? 0 : wig(i - 1, j)) + wig(i, j) * f1 * (f2 + cosBeta[j]);
        }
    }
}

PSOFFT_NAMESPACE_END

#endif /* fn_dwt.hpp */
