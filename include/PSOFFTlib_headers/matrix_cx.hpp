//
//  matrix_cx.hpp
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

#ifndef PSOFFTlib_matrix_cx_hpp
#define PSOFFTlib_matrix_cx_hpp

PSOFFT_BEGIN

/*!
 * @brief       Matrix specialization for complex data.
 * @details     A matrix that contains only complex values on each
 *              entry.  The matrix class uses internally an dynamic allocated array
 *              to store values. The memory is organized in  column-major
 *              order to provide interchangability with the BLAS and LAPACK
 *              Fortran libraries.
 *
 * @tparam      T An element type which represents a number that provides all common
 *              mathmatical operations.
 *
 * @since       0.0.1
 * 
 * @todo        Implement all other additional operators for 'lhs * rhs' outside
 *              of class definition
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        13.05.15
 *
 * @ingroup     matrix
 */
template< typename T >
class
matrix< complex< T >, if_pod_type< T > >
{
public:
    // type declarations
    typedef T pod_type;             //!< The POD type of matrix elements
    typedef typename smart_array< complex< pod_type > >::iterator lin_iterator;
    
    // ivars
    const smart_array< complex< pod_type > > mem;   //!< Matrix memory
    const size_t rows;                              //!< Number of matrix rows
    const size_t cols;                              //!< Number of matrix columns
    
    // methods
    inline                                     ~matrix();
    inline                                      matrix();
    
    inline                                      matrix(const size_t& m, const size_t& n);
    inline                                      matrix(const size_t& mn);
    inline                                      matrix(const size_t& m, const size_t& n, const complex< pod_type >& initial);
    inline                                      matrix(const matrix< complex< pod_type > >& A);
    inline                                      matrix(const matrix< pod_type >& A);
    
    inline const matrix< complex< pod_type > >& operator=(const matrix< pod_type >& A);
    inline const matrix< complex< pod_type > >& operator=(const matrix< complex< pod_type > >& A);
    
    inline const matrix< complex< pod_type > >& operator*=(const vector< pod_type >& v);
    inline const matrix< complex< pod_type > >& operator*=(const vector< complex< pod_type > >& v);
    
    inline       matrix< complex< pod_type > >  operator+();
    inline       matrix< complex< pod_type > >  operator-();
    
    inline       matrix< complex< pod_type > >  operator*(const pod_type& rhs);
    inline       matrix< complex< pod_type > >  operator*(const complex< pod_type >& rhs);
    
    inline       vector< complex< pod_type > >  operator*(const vector< pod_type >& v);
    inline       vector< complex< pod_type > >  operator*(const vector< complex< pod_type > >& v);
    
    inline       matrix< complex< pod_type > >& operator*=(const pod_type& rhs);
    inline       matrix< complex< pod_type > >& operator*=(const complex< pod_type >& rhs);
    
    inline       complex< pod_type >&           operator()(const size_t& i, const size_t& j);
    inline const complex< pod_type >&           operator()(const size_t& i, const size_t& j) const;
    
    inline       void                           transpose();
    
    inline       lin_iterator                   begin();
    inline       lin_iterator                   end();
};

/*!
 * @brief           The deconstructor to delete a complex matrix and its contents properly from
 *                  RAM.
 * @details         Deletes the complex matrix and frees its allocated dynamic mamemory.
 */
template< typename T >
inline
matrix< complex< T >, if_pod_type< T > >::~matrix()
{}

/*!
 * @brief           Default constructor to build a complex matrix with no rows and columns.
 * @details         Constructs a complex matrix that has no columns and rows and is declared
 *                  as not initialized.
 */
template< typename T >
inline
matrix< complex< T >, if_pod_type< T > >::matrix()
    : rows(0)
    , cols(0)
{}

/*!
 * @brief           Constructs a complex matrix with given parameters
 * @details         Constructs the complex matrix by allocating enough memory to store \f$M\times N\f$
 *                  values. The content of the complex matrix is initialized with \f$0 + 0i\f$
 *
 * @param[in]       m Number of rows in the constructed complex matrix
 * @param[in]       n Number of columns in the constructed complex matrix
 */
template< typename T >
inline
matrix< complex< T >, if_pod_type< T > >::matrix(const size_t& m, const size_t& n)
    : rows(m)
    , cols(n)
{
    access::rw(mem) = smart_array< complex< pod_type > >(rows * cols);
}

/*!
 * @brief           Constructs a square complex matrix with given dimension \f$N\times N\f$
 * @details         Constructs the complex matrix by allocating enough memory to store \f$N\times N\f$
 *                  values. The content of the complex matrix is initialized with \f$0+0i\f$.
 *
 * @param[in]       mn Number of rows and columns in the constructed complex matrix
 */
template< typename T >
inline
matrix< complex< T >, if_pod_type< T > >::matrix(const size_t& mn)
    : rows(mn)
    , cols(mn)
{
    access::rw(mem) = smart_array< complex< pod_type > >(rows * cols);
}

/*!
 * @brief           Constructs a complex matrix with given dimensions and given intial element value
 * @details         Constructs the complex matrix by allocating enough memory to store \f$M\times N\f$
 *                  values. Each value in the complex matrix will be initialized with given initial value.
 *
 * @param[in]       m Number of rows in the constructed complex matrix
 * @param[in]       n Number of columns in the constructed complex matrix
 * @param[in]       initial Initial value for each complex matrix element
 */
template< typename T >
inline
matrix< complex< T >, if_pod_type< T > >::matrix(const size_t& m, const size_t& n, const complex< pod_type >& initial)
    : rows(m)
    , cols(n)
{
    access::rw(mem) = smart_array< complex< pod_type > >(rows * cols);
    
    if (rows * cols > 0)
    {
        for (int i = 0; i < rows * cols; ++i)
        {
            access::rw(mem[i]) = initial;
        }
    }
}

/*!
 * @brief           A copy constructor for copying a given complex matrix
 * @details         Copies the contents of the given complex matrix to build a new complex matrix
 *                  with the same state.
 *
 * @param[in]       A The complex matrix that is supposed to be copied.
 */
template< typename T >
inline
matrix< complex< T >, if_pod_type< T > >::matrix(const matrix< complex< pod_type > >& A)
    : rows(A.rows)
    , cols(A.cols)
{
    access::rw(mem) = A.mem;
}

/*!
 * @brief           A copy constructor for copying a given matrix
 * @details         Copies the contents of the given matrix to build a new complex matrix
 *                  with the same state.
 *
 * @param[in]       A The matrix that is supposed to be copied.
 */
template< typename T >
inline
matrix< complex< T >, if_pod_type< T > >::matrix(const matrix< pod_type >& A)
    : rows(A.rows)
    , cols(A.cols)
{
    access::rw(mem) = smart_array< complex< pod_type > >(rows * cols);
    
    for (int i = 0; i < rows * cols; ++i)
    {
        access::rw(mem[i]) = complex< pod_type >(A.mem[i], 0);
    }
}

/*!
 * @brief           The assignment operator to copy matrices.
 * @details         Copies the given matrix \f$A\f$ into a new complex matrix
 *
 * @param[in]       A The matrix that is supposed to be copied
 *
 * @return          A new matrix that has the same state than \f$A\f$ and containing
 *                  the same values.
 */
template< typename T >
inline
const matrix< complex< T > >& matrix< complex< T >, if_pod_type< T > >::operator=(const matrix< pod_type >& A)
{
    access::rw(rows) = A.rows;
    access::rw(cols) = A.cols;
    access::rw(mem)  = smart_array< complex< pod_type > >(rows * cols);
    
    if (rows * cols > 0)
    {
        for (int i = 0; i < rows * cols; ++i)
        {
            access::rw(mem[i]) = complex< pod_type >(A.mem[i], 0);
        }
    }
    
    return *this;
}

/*!
 * @brief           The assignment operator to copy a complex matrices.
 * @details         Copies the given complex matrix \f$A\f$ into a new complex matrix
 *
 * @param[in]       A The complex matrix that is supposed to be copied
 *
 * @return          A new complex matrix that has the same state than \f$A\f$ and containing
 *                  the same values.
 */
template< typename T >
inline
const matrix< complex< T > >& matrix< complex< T >, if_pod_type< T > >::operator=(const matrix< complex< pod_type > >& A)
{
    if ( this == &A )
    {
        return *this;
    }
    
    access::rw(rows) = A.rows;
    access::rw(cols) = A.cols;
    access::rw(mem)  = A.mem;
    
    return *this;
}

/*!
 * @brief           The multiplication assignment operator for matrix-complex vector product.
 * @details         Multiplying the complex vector \f$v\f$ to the current matrix and stores the
 *                  result in the current matrix.
 *
 * @param[in]       v The complex vector that is supposed to be multiplied to the current matrix
 *
 * @return          The reference to the current matrix that contains the result.
 */
template< typename T >
inline
const matrix< complex< T > >& matrix< complex< T >, if_pod_type< T > >::operator*=(const vector< complex< pod_type > >& v)
{
    typedef vector< complex< pod_type > > cx_T_vector;
    psofft_error(v.type == cx_T_vector::ROW || cols != v.size, "%s", "Dimension mismatch in complex matrix-vector multiplication.");
    
    // create new memory array
    smart_array< complex< pod_type > > new_mem(rows);
    
    // set each value to 0
    for (typename smart_array< complex< pod_type > >::iterator e = new_mem.begin(); e != new_mem.end(); ++e)
    {
        *e = 0;
    }
    
    // do multiplication
    size_t i, j;
    for (j = 0; j < cols; ++j)
    {
        for (i = 0; i < rows; ++i)
        {
            new_mem[i] +=  access::rw(mem[j * rows + i]) * v[j];
        }
    }
    
    // adjust size and memory
    access::rw(cols) = 1;
    access::rw(mem)  = new_mem;
    
    // return reference to the current matrix
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
matrix< complex< T > > matrix< complex< T >, if_pod_type< T > >::operator+()
{
    matrix< complex< pod_type > > C(rows, cols);
    access::rw(C.mem) = mem;
    
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
matrix< complex< T > > matrix< complex< T >, if_pod_type< T > >::operator-()
{
    matrix< complex< pod_type > > C(rows, cols);
    
    size_t i;
    for (i = 0; i < rows * cols; ++i)
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
inline matrix< complex< T > > matrix< complex< T >, if_pod_type< T > >::operator*(const pod_type& rhs)
{
    matrix< complex< pod_type > > C(rows, cols);
    
    size_t i;
    for (i = 0; i < rows * cols; ++i)
    {
        access::rw(C.mem[i]) = mem[i] * complex< pod_type >(rhs, 0);
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
matrix< complex< T > > matrix< complex< T >, if_pod_type< T > >::operator*(const complex< pod_type >& rhs)
{
    matrix< complex< pod_type > > C(rows, cols);
    
    size_t i;
    for (i = 0; i < rows * cols; ++i)
    {
        access::rw(C.mem[i]) = mem[i] * rhs;
    }
    
    return C;
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
matrix< complex< T > >& matrix< complex< T >, if_pod_type< T > >::operator*=(const pod_type& rhs)
{
    size_t i;
    for (i = 0; i < rows * cols; ++i)
    {
        access::rw(mem[i]) *= complex< pod_type >(rhs, 0);
    }
    return *this;
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
matrix< complex< T > >& matrix< complex< T >, if_pod_type< T > >::operator*=(const complex< pod_type >& rhs)
{
    size_t i;
    for (i = 0; i < rows * cols; ++i)
    {
        access::rw(mem[i]) *= rhs;
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
 * @return          The pointer to matrix element in row \f$i\f$ and column \f$j\f$ (\f$a_{ij}\f$).
 */
template< typename T >
inline
complex< T >& matrix< complex< T >, if_pod_type< T > >::operator()(const size_t& i, const size_t& j)
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
const complex< T >& matrix< complex< T >, if_pod_type< T > >::operator()(const size_t& i, const size_t& j) const
{
    return mem[j * rows + i];
}

/*!
 * @brief           Transposing the current matrix.
 * @details         Swapping the elements of the current matrix according to
 *                  \f[
 *                      A_{ij} = A_{ji}
 *                  \f]
 * @todo            Implement in-place transposing
 */
template< typename T >
inline
void matrix< complex< T >, if_pod_type< T > >::transpose()
{
    smart_array< complex< pod_type > > tmp_mem(rows * cols);
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
    
    access::rw(rows) = tmp_r;
    access::rw(cols) = tmp_c;
    access::rw(mem)  = tmp_mem;
}

/*!
 * @brief           Outstream operator overload for complex matrices.
 * @details         The out-steam operator is used to print the matrix
 *                  in a nice form over the std::cout stream.
 * @param[in,out]   o The stream object
 * @param[in]       A The matrix that should be printed
 *
 * @return          The reference to the given out-stream.
 *
 * @ingroup         matrix
 */
template< typename T >
std::ostream& operator<<(std::ostream& o, const matrix< complex< T > >& A)
{
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl;
    
    int width   = 20;
    auto format = std::fixed;
    
    if ( different_type< T, float >::value && different_type< T, double >::value && different_type< T, long double >::value )
    {
        width = 10;
    }
    
    // check values
    size_t i, j;
    for (i = 0; i < A.rows; ++i)
    {
        for (j = 0; j < A.cols; ++j)
        {
            complex< T > c = A(i, j);
            if (std::abs(c.re) >= 10 || std::abs(c.im) >= 10)
            {
                width   = 22;
                format  = std::fixed;
                
                if ( different_type< T, float >::value && different_type< T, double >::value && different_type< T, long double >::value )
                {
                    width = 12;
                }
            }
            
            if (std::abs(c.re) >= 100 || std::abs(c.im) >= 100)
            {
                width   = 24;
                format  = std::fixed;
                
                if ( different_type< T, float >::value && different_type< T, double >::value && different_type< T, long double >::value )
                {
                    width = 14;
                }
            }
            
            if (std::abs(c.re) >= 1000 || std::abs(c.im) >= 1000)
            {
                width   = 28;
                format  = std::scientific;
                
                if ( different_type< T, float >::value && different_type< T, double >::value && different_type< T, long double >::value )
                {
                    width = 18;
                }
            }
        }
    }
    
    // prepare output and print
    for (i = 0; i < A.rows; ++i)
    {
        for (j = 0; j < A.cols; ++j)
        {
            // get entry
            complex< T > c = A(i, j);
            
            // create string
            std::ostringstream val;
            
            // add real value to string
            val << format << std::setprecision(4) << c.re;
            val << (c.im < 0 ? " - " : " + ") << (c.im == 0 ?  0 : std::abs(c.im)) << "i";
            
            // get string from stream
            std::string str = val.str();
            
            // set filling character
            o << std::setfill(' ') << std::right << std::setw(width) << str;
        }
        
        // line break for next line
        o << std::endl;
    }
    
    std::cout.flags( f );
    return o;
}

PSOFFT_END

#endif /* matrix_cx.hpp */
