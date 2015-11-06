//
//  fn_flip.hpp
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

#ifndef PSOFFTlib_fn_flip_hpp
#define PSOFFTlib_fn_flip_hpp

PSOFFT_BEGIN

/*!
 * @brief           Flips a given matrix by mirroring the columns.
 * @details         Flips a matrix so that 
 *                  \f$M = (v_0,v_1,v_2,\dots,v_{n-2}, v_{n-1}, v_{n})\f$
 *                  becomes
 *                  \f$M = (v_n,v_{n-1},v_{n-2},\dots,v_2,v_1,v_0)\f$
 *
 * @param[in,out]   mat The matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            16.06.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename T >
inline
void_number_type< T > fliplr(matrix< T >& mat)
{
    // define indices
    size_t j, k;
    
    // iterate over half of columns
    for (j = 0; j < mat.cols / 2; ++j)
    {
        // iterate over rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (k = 0; k < mat.rows; ++k)
        {
            // temporary element copy
            T tmp = mat(k, j);
            
            // swap elements
            mat(k, j)                = mat(k, mat.cols - j - 1);
            mat(k, mat.cols - j - 1) = tmp;
        }
    }
}

/*!
 * @brief           Flips a given matrix by inverting the column order and negating
 *                  each element in every second even row.
 * @details         Flips a matrix so that
 *                  \f$M = (v_0,v_1,v_2,\dots,v_{n-2}, v_{n-1}, v_{n})\f$
 *                  becomes
 *                  \f$M = (v_n,v_{n-1},v_{n-2},\dots,v_2,v_1,v_0)\f$
 *                  with negated elements in every even row.
 *
 * @param[in,out]   mat The matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            16.06.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename T >
inline
void_number_type< T > fliplr_ne2nderow(matrix< T >& mat)
{
    // define indices
    size_t j, k;
    
    // iterate over half of columns
    for (j = 0; j < mat.cols / 2; ++j)
    {
        // iterate over rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (k = 0; k < mat.rows; ++k)
        {
            // temporary element copy
            T tmp = mat(k, j);
            
            // swap elements
            if (k & 1)
            {
                mat(k, j)                = mat(k, mat.cols - j - 1);
                mat(k, mat.cols - j - 1) = tmp;
            }
            else
            {
                mat(k, j)                = -mat(k, mat.cols - j - 1);
                mat(k, mat.cols - j - 1) = -tmp;
            }
        }
    }
}

/*!
 * @brief           Flips a given matrix by inverting the column order and negating
 *                  each element in every second odd row.
 * @details         Flips a matrix so that
 *                  \f$M = (v_0,v_1,v_2,\dots,v_{n-2}, v_{n-1}, v_{n})\f$
 *                  becomes
 *                  \f$M = (v_n,v_{n-1},v_{n-2},\dots,v_2,v_1,v_0)\f$
 *                  with negated element in every odd row.
 *
 * @param[in,out]   mat The matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            16.06.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename T >
inline
void_number_type< T > fliplr_ne2ndorow(matrix< T >& mat)
{
    // define indices
    size_t j, k;
    
    // iterate over half of columns
    for (j = 0; j < mat.cols / 2; ++j)
    {
        // iterate over rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (k = 0; k < mat.rows; ++k)
        {
            // temporary element copy
            T tmp = mat(k, j);
            
            // swap elements
            if (k & 1)
            {
                mat(k, j)                = -mat(k, mat.cols - j - 1);
                mat(k, mat.cols - j - 1) = -tmp;
            }
            else
            {
                mat(k, j)                = mat(k, mat.cols - j - 1);
                mat(k, mat.cols - j - 1) = tmp;
            }
        }
    }
}

/*!
 * @brief           Flips a given matrix inverting the row order.
 * @details         Flips a matrix so that
 *                  \f$M = \begin{pmatrix}v_0\\ v_1\\ v_2\\ \dots\\ v_{n-2}\\ v_{n-1}\\ v_{n}\end{pmatrix}\f$
 *                  becomes
 *                  \f$M = \begin{pmatrix}v_n\\ v_{n-1}\\ v_{n-2}\\ \dots\\ v_2\\ v_1\\ v_0\end{pmatrix}\f$
 *
 * @param[in,out]   mat The matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            03.07.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename T >
inline
void_number_type< T > flipud(matrix< T >& mat)
{
    // define indices
    size_t j, k;
    
    // iterate over cols
    for (k = 0; k < mat.cols; ++k)
    {
        // iterate over half of rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (j = 0; j < mat.rows / 2; ++j)
        {
            // temporary element copy
            T tmp = mat(j, k);
            // swap elements
            mat(j, k)                = mat(mat.rows - j - 1, k);
            mat(mat.rows - j - 1, k) = tmp;
        }
    }
}

/*!
 * @brief           Flips a given matrix by inverting the row order and negating
 *                  each element in every second even column.
 * @details         Flips a matrix so that
 *                  \f$M = (v_0,v_1,v_2,\dots,v_{n-2}, v_{n-1}, v_{n})\f$
 *                  becomes
 *                  \f$M = (v_n,v_{n-1},v_{n-2},\dots,v_2,v_1,v_0)\f$
 *                  with negated elements in every second even column.
 *
 * @param[in,out]   mat The matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            16.06.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename T >
inline
void_number_type< T > flipud_ne2ndecol(matrix< T >& mat)
{
    // define indices
    size_t j, k;
    
    // iterate over cols
    for (k = 0; k < mat.cols; ++k)
    {
        // iterate over half of rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (j = 0; j < mat.rows / 2; ++j)
        {
            // temporary element copy
            T tmp = mat(j, k);
            
            // swap elements
            if (k & 1)
            {
                mat(j, k)                = mat(mat.rows - j - 1, k);
                mat(mat.rows - j - 1, k) = tmp;
            }
            else
            {
                mat(j, k)                = -mat(mat.rows - j - 1, k);
                mat(mat.rows - j - 1, k) = -tmp;
            }
        }
    }
}

/*!
 * @brief           Flips a given matrix by inverting the row order and negating
 *                  each element in every second odd column.
 * @details         Flips a matrix so that
 *                  \f$M = (v_0,v_1,v_2,\dots,v_{n-2}, v_{n-1}, v_{n})\f$
 *                  becomes
 *                  \f$M = (v_n,v_{n-1},v_{n-2},\dots,v_2,v_1,v_0)\f$
 *                  with negated elements in every second odd column.
 *
 * @param[in,out]   mat The matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            16.06.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename T >
inline
void_number_type< T > flipud_ne2ndocol(matrix< T >& mat)
{
    // define indices
    size_t j, k;
    
    // iterate over cols
    for (k = 0; k < mat.cols; ++k)
    {
        // iterate over half of rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (j = 0; j < mat.rows / 2; ++j)
        {
            // temporary element copy
            T tmp = mat(j, k);
            
            // swap elements
            if (k & 1)
            {
                mat(j, k)                = -mat(mat.rows - j - 1, k);
                mat(mat.rows - j - 1, k) = -tmp;
            }
            else
            {
                mat(j, k)                = mat(mat.rows - j - 1, k);
                mat(mat.rows - j - 1, k) = tmp;
            }
        }
    }
}

PSOFFT_END

#endif /* fn_flip.hpp */
