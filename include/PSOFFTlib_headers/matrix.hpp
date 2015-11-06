//
//  matrix.hpp
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

#ifndef PSOFFTlib_matrix_hpp
#define PSOFFTlib_matrix_hpp

PSOFFT_BEGIN

/*!
 * @brief       Collection of classes and functions for matrices for mathematical
 *              purposes.
 * @defgroup    matrix Matrix
 * @{
 */

/*!
 * @brief       Matrix class with common operations.
 * @details     The matrix class uses internally an dynamic allocated array
 *              to store values. The memory is organized in  column-major
 *              order to provide interchangability with the BLAS and LAPACK
 *              Fortran libraries.
 *
 * @tparam      T An element type which represents a number that provides all common
 *              mathmatical operations.
 *
 * @since       0.0.1
 *
 * @todo        Implement matrix transposition with flags to change setting
 *              and getting behavior of the ()-operator. This will save some
 *              extra computation
 * @todo        Implement injection functionality and token for the <<-Operator.
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        12.01.15
 */
template< typename T >
class
matrix< T, if_pod_type< T > >
{
public:
    typedef T pod_type;     //!< The POD type of matrix elements
    
    // ivars
    const size_t rows;      //!< Number of matrix rows
    const size_t cols;      //!< Number of matrix columns
    
    const pod_type* mem;    //!< Matrix memory
    
    // methods
    inline                                     ~matrix();
    inline                                      matrix();
    
    inline                                      matrix(const size_t& m, const size_t& n);
    inline                                      matrix(const size_t& mn);
    inline                                      matrix(const size_t& m, const size_t& n, const pod_type& initial);
    
    inline                                      matrix(const matrix< pod_type >& A);
    inline                                      matrix(matrix< pod_type >&& A);
    
    inline       vector< pod_type >             operator*(const vector< pod_type >& v);
    inline       vector< complex< pod_type > >  operator*(const vector< complex< pod_type > >& v);
    
    inline const matrix< pod_type >&            operator=(const matrix< pod_type >& A);
    inline const matrix< pod_type >&            operator=(matrix< pod_type >&& A);
    
    inline const matrix< pod_type >&            operator*=(const vector< pod_type >& v);
    inline       matrix< pod_type >&            operator*=(const pod_type& rhs);
    
    inline       matrix< pod_type >             operator+();
    inline       matrix< pod_type >             operator-();
    
    inline       matrix< pod_type >                             operator*(const pod_type& rhs);
    inline       matrix< complex< pod_type > >  operator*(const complex< pod_type >& rhs);
    
    inline       pod_type&                      operator()(const size_t& i, const size_t& j);
    inline const pod_type&                      operator()(const size_t& i, const size_t& j) const;
    
    inline       void                           transpose();
};

/*!
 * @brief           Default constructor to build a matrix with no rows and columns.
 * @details         Constructs a matrix that has no columns and rows and is declared
 *                  as not initialized.
 */
template< typename T >
inline
matrix< T, if_pod_type< T > >::matrix()
    : rows(0)
    , cols(0)
    , mem(nullptr)
{}

/*!
 * @brief           The deconstructor to delete a matrix and its contents properly from
 *                  RAM.
 * @details         Deletes the matrix and frees its allocated dynamic mamemory.
 */
template< typename T >
inline
matrix< T, if_pod_type< T > >::~matrix()
{
    delete [] mem;
}

/*!
 * @brief           Constructs a matrix with given parameters
 * @details         Constructs the matrix by allocating enough memory to store \f$M\times N\f$
 *                  values. The content of the matrix is not initialized and can sometimes contain
 *                  random values that were written to the allocated memory before.
 *
 * @param[in]       m Number of rows in the constructed matrix
 * @param[in]       n Number of columns in the constructed matrix
 */
template< typename T >
inline
matrix< T, if_pod_type< T > >::matrix(const size_t& m, const size_t& n)
    : rows(m)
    , cols(n)
{
    size_t cap  = m * n;
    mem         = new T[cap];
}

/*!
 * @brief           Constructs a square matrix with given dimension \f$N\times N\f$
 * @details         Constructs the matrix by allocating enough memory to store \f$N\times N\f$
 *                  values. The content of the matrix is not initialized and can sometimes contain
 *                  random values that were written to the allocated memory before.
 *
 * @param[in]       mn Number of rows and columns in the constructed matrix
 */
template< typename T >
inline
matrix< T, if_pod_type< T > >::matrix(const size_t& mn)
    : rows(mn)
    , cols(mn)
{
    size_t cap  = mn * mn;
    mem         = new T[cap];
}

/*!
 * @brief           Constructs a matrix with given dimensions and given intial element value
 * @details         Constructs the matrix by allocating enough memory to store \f$M\times N\f$
 *                  values. Each value in the matrix will be initialized with given initial value.
 *
 * @param[in]       m Number of rows in the constructed matrix
 * @param[in]       n Number of columns in the constructed matrix
 * @param[in]       initial Initial value for each matrix element
 */
template< typename T >
inline
matrix< T, if_pod_type< T > >::matrix(const size_t& m, const size_t& n, const T& initial)
    : rows(m)
    , cols(n)
{
    size_t cap  = m * n;
    mem         = new T[cap];
    
    if (cap > 0)
    {
        if (initial == 0 || initial == -1)
        {
            // fastest possible assembler routine
            memset(access::rwp(mem), initial, cap * sizeof(T));
        }
        else
        {
            std::fill(access::rwp(mem), access::rwp(mem) + cap, initial);
        }
    }
}

/*!
 * @brief           A copy constructor for copying a given matrix
 * @details         Copies the contents of the given matrix to build a new matrix
 *                  with the same state.
 *
 * @param[in]       A The matrix that is supposed to be copied.
 */
template< typename T >
inline
matrix< T, if_pod_type< T > >::matrix(const matrix< T >& A)
    : rows(A.rows)
    , cols(A.cols)
{
    size_t cap  = rows * cols;
    mem         = new T[cap];
    
    if (cap > 0)
    {
        memcpy(access::rwp(mem), access::rwp(A.mem), cap * sizeof(T));
    }
}

/*!
 * @brief           A move constructor to copy an r-value matrix into a new matrix
 * @details         Takes the memory of the r-value matrix \f$A\f$ and copies all other
 *                  instance variables. The resulting matrix is in the same state
 *                  than the original matrix. A move constructor has better performance
 *                  by assigning matrices that were calculated in one expression.
 *
 * @param[in,out]   A The r-value matrix \f$A\f$ which content should be moved to a
 *                  new matrix.
 */
template< typename T >
inline
matrix< T, if_pod_type< T > >::matrix(matrix< T >&& A)
    : rows(A.rows)
    , cols(A.cols)
{
    const T* tmp = mem;
    mem           = A.mem;
    A.mem         = tmp;
}

/*!
 * @brief           The assignment operator to copy matrices.
 * @details         Copies the given matrix \f$A\f$ into a new matrix
 *
 * @param[in]       A The matrix that is supposed to be copied
 *
 * @return          A new matrix that has the same state than \f$A\f$ and containing
 *                  the same values.
 */
template< typename T >
inline
const matrix< T >&  matrix< T, if_pod_type< T > >::operator=(const matrix< T >& A)
{
    if ( this == &A )
    {
        return *this;
    }
    
    rows        = A.rows;
    cols        = A.cols;
    size_t cap  = rows * cols;
    
    delete [] mem;
    
    mem = new T[cap];
    
    if (cap > 0)
    {
        memcpy(mem, A.mem, cap * sizeof(T));
    }
    
    return *this;
}

/*!
 * @brief           The assignment operator to copy matrices by moving its contents.
 * @details         Moving the contents of a r-value matrix into a new matrix
 *
 * @param[in]       A The r-value matrix that is supposed to be copied
 *
 * @return          A new matrix that has the same state than \f$A\f$ and containing
 *                  the same values.
 */
template< typename T >
inline
const matrix< T >&  matrix< T, if_pod_type< T > >::operator=(matrix< T >&& A)
{
    if ( this == &A )
    {
        return *this;
    }
    
    access::rw(rows) = A.rows;
    access::rw(cols) = A.cols;
    
    const T* tmp = mem;
    mem           = A.mem;
    A.mem         = tmp;
    
    return *this;
}

/*!
 * @brief           The plus sign operator.
 * @details         Copies the matrix applying
 *                  \f{eqnarray*}{
 *                      B_{ij} = +A_{ij}
 *                  \f}
 *
 * @return          A matrix containing a copy of the current matrix.
 */
template< typename T >
inline
matrix< T > matrix< T, if_pod_type< T > >::operator+()
{
    matrix< T > C(rows, cols);
    memcpy(C.mem, mem, rows * cols * sizeof(T));
    
    return C;
}

/*!
 * @brief           The minus sign operator.
 * @details         Copies the matrix applying
 *                  \f{eqnarray*}{
 *                      B_{ij} = -A_{ij}
 *                  \f}
 *
 * @return          A matrix containing the negative values of the current matrix.
 */
template< typename T >
inline
matrix< T > matrix< T, if_pod_type< T > >::operator-()
{
    matrix< T > C(rows, cols);
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        access::rw(C.mem[i]) = - mem[i];
    }
    
    return C;
}

/*!
 * @brief           The multiply operator for scaling the matrix elements.
 * @details         The scalar will be multiplied with each element in the current
 *                  matrix.
 *
 * @param[in]       rhs The scalar value
 *
 * @return          A new matrix containing the result of multiplication
 *                  \f{eqnarray*}{
 *                      \lambda \left(\begin{array}{c c c c}
 *                          a_{11} & a_{12} & \cdots & a_{1n}\\
 *                          a_{21} & a_{22} & \cdots & a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} &  a_{m2} & \cdots & a_{mn}
 *                      \end{array}\right) =
 *                      \left(\begin{array}{c c c c}
 *                          \lambda a_{11} & \lambda a_{12} & \cdots & \lambda a_{1n}\\
 *                          \lambda a_{21} & \lambda a_{22} & \cdots & \lambda a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          \lambda a_{m1} & \lambda a_{m2} & \cdots & \lambda a_{mn}
 *                      \end{array}\right)
 *                  \f}
 */
template< typename T >
inline
matrix< T > matrix< T, if_pod_type< T > >::operator*(const T& rhs)
{
    // result matrix
    matrix< T > C(rows, cols);
    
    // capacity
    size_t cap = rows * cols;
    
    // iterate over temporary and internal memory and scale data
    for (const T *e1 = mem, *e2 = C.mem; e1 != mem + cap; ++e1, ++e2)
    {
        access::rw(*e2) = *e1 * rhs;
    }
    
    // return result
    return C;
}

/*!
 * @brief           The multiply operator for scaling the matrix elements.
 * @details         The scalar will be multiplied with each element in the current
 *                  matrix.
 *
 * @param[in]       rhs The scalar value
 *
 * @return          A new matrix containing the result of multiplication
 *                  \f{eqnarray*}{
 *                      \lambda \left(\begin{array}{c c c c}
 *                          a_{11} & a_{12} & \cdots & a_{1n}\\
 *                          a_{21} & a_{22} & \cdots & a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} &  a_{m2} & \cdots & a_{mn}
 *                      \end{array}\right) =
 *                      \left(\begin{array}{c c c c}
 *                          \lambda a_{11} & \lambda a_{12} & \cdots & \lambda a_{1n}\\
 *                          \lambda a_{21} & \lambda a_{22} & \cdots & \lambda a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          \lambda a_{m1} & \lambda a_{m2} & \cdots & \lambda a_{mn}
 *                      \end{array}\right)
 *                  \f}
 */
template< typename T >
inline
matrix< complex< T > > matrix< T, if_pod_type< T > >::operator*(const complex< T >& rhs)
{
    matrix< complex< T > > C(rows, cols);
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        C.mem = complex< T >(mem[i], 0) * rhs;
    }
    
    return C;
}

/*!
 * @brief           The matrix vector multiplication operator.
 * @details         Multiplies the current matrix with the given vector.
 *
 * @param[in]       v The vector which is supposed to be multiplied with the current
 *                  matrix
 *
 * @return          The result of the mutliplication
 */
template< typename T >
inline
vector< complex< T > > matrix< T, if_pod_type< T > >::operator*(const vector< complex< T > >& v)
{
    if (cols != v.size || v.type == vector< complex< T > >::type::ROW)
    {
        psofft_error("%s", "dimension mismatch in matrix-complex vector multiplication.");
    }
    
    vector< complex< T > > result(rows, 0, v.type);
    
    size_t i, j;
    for (i = 0; i < cols; ++i)
    {
        for (j = 0; j < rows; ++j)
        {
            result[j] += complex< T >(mem[i * rows + j], 0) * v[i];
        }
    }
    
    return result;
}

/*!
 * @brief           The multiply assignment operator for a scalar value.
 * @details         Multiplies the scalar value to each entry of the current matrix.
 *
 * @param[in]       rhs The scalar value that is supposed to be multiplied to each entry.
 *
 * @return          The reference to the current matrix.
 */
template< typename T >
inline
matrix< T >& matrix< T, if_pod_type< T > >::operator*=(const T& rhs)
{
    size_t cap = rows * cols;
    for (const T* e = mem; e != mem + cap; ++e)
    {
        access::rw(*e) *= rhs;
    }
    
    return *this;
}

/*!
 * @brief           Indexing operator to access the specif matrix element.
 * @details         The indexing operator can be used to set a specific element of the
 *                  current matrix by using two indices \f$m\f$ and \f$n\f$.
 *
 * @param[in]       i The row of the needed element
 * @param[in]       j The column of the needed element
 *
 * @return      The pointer to matrix element in row \f$i\f$ and column \f$j\f$ (\f$a_{ij}\f$).
 */
template< typename T >
inline
T& matrix< T, if_pod_type< T > >::operator()(const size_t& i, const size_t& j)
{
    return access::rw(mem[j * rows + i]);
}

/*!
 * @brief           Indexing operator to get the specific matrix element.
 * @details         The indexing operator can be used to get a specific element of the
 *                  current matrix by using two indices \f$m\f$ and \f$n\f$.
 *
 * @param[in]       i The row of the needed element
 * @param[in]       j The column of the needed element
 *
 * @return          The matrix element in row \f$i\f$ and column \f$j\f$ (\f$a_{ij}\f$).
 */
template< typename T >
inline
const T& matrix< T, if_pod_type< T > >::operator()(const size_t& i, const size_t& j) const
{
    return mem[j * rows + i];
}

/*!
 * @brief           Transposing the current matrix.
 * @details         Swapping the elements of the current matrix according to
 *                  \f[
 *                      A_{ij} = A_{ji}
 *                  \f]
 * @todo            Implement in-place transpose
 */
template< typename T >
inline
void matrix< T, if_pod_type< T > >::transpose()
{
    T *tmp_mem  = new T[rows * cols];
    size_t tmp_r = cols;
    size_t tmp_c = rows;
    
    size_t i, j;
    for (i = 0; i < rows; ++i)
    {
        for (j = 0; j < cols; ++j)
        {
            tmp_mem[i * cols + j] = mem[j * rows + i];
        }
    }
    
    delete [] mem;
    access::rw(rows) = tmp_r;
    access::rw(cols) = tmp_c;
    mem              = tmp_mem;
}

/*!
 * @brief           Outstream operator overload.
 * @details         The out-steam operator is used to print the matrix
 *                  in a nice form over the std::cout stream.
 * @param[in,out]   o The stream object
 * @param[in]       A The matrix that should be printed
 *
 * @return          The reference to the given out-stream.
 */
template< typename S >
std::ostream& operator<<(std::ostream& o, const matrix< S >& A)
{
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl;
    
    int width   = 10;
    auto format = std::fixed;
    
    // reduce size for integers
    if ( different_type< S, float >::value == false && different_type< S, double >::value && different_type< S, long double >::value )
    {
        width = 5;
    }
    
    // check values
    size_t i, j;
    for (i = 0; i < A.rows; ++i)
    {
        for (j = 0; j < A.cols; ++j)
        {
            S r = A(i, j);
            if (std::abs(r) >= 10)
            {
                width   = 11;
                format  = std::fixed;
                
                if ( different_type< S, float >::value == false && different_type< S, double >::value && different_type< S, long double >::value )
                {
                    width = 6;
                }
            }
            
            if (std::abs(r) >= 100)
            {
                width   = 12;
                format  = std::fixed;
                
                if ( different_type< S, float >::value == false && different_type< S, double >::value && different_type< S, long double >::value )
                {
                    width = 7;
                }
            }
            
            if (std::abs(r) >= 1000)
            {
                width   = 14;
                format  = std::scientific;
                
                if ( different_type< S, float >::value == false && different_type< S, double >::value && different_type< S, long double >::value )
                {
                    width = 10;
                }
            }
        }
    }
    
    // setting decimal precesion
    for (i = 0; i < A.rows; ++i)
    {
        for (j = 0; j < A.cols; ++j)
        {
            // get entry
            S val = A(i, j);
            
            // create string
            o << std::setw(width);      // setting fixed width for the number
            o << std::setprecision(4);  // setting number precision
            o << std::setfill(' ');     // fill space with white spaces
            o << std::right;            // setting right alignment of number
            o << format;                // setting correct number formatting
            o << val;                   // print value
        }
        o << std::endl;
    }
    
    std::cout.flags( f );
    return o;
}

/*!
 * @}
 */

PSOFFT_END

#endif /* matrix.hpp */
