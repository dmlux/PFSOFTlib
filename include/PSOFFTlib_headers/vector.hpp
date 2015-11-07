//
//  vector.hpp
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

#ifndef PSOFFTlib_vector_hpp
#define PSOFFTlib_vector_hpp

PSOFFT_BEGIN

/*!
 * @brief       Collection of classes and functions for vectors for mathematical 
 *              purposes.
 * @defgroup    vector Vector
 * @{
 */
 
/*!
 * @brief       Vector class with special operations and methods.
 * @details     The vector class uses internally an dynamic allocated array to store 
 *              values. The memory is organized in as linear memory. It provides a 
 *              set of methods and operators to provide a useful vector class for 
 *              mathematical purposes.
 *
 *
 * @tparam      T An element type which represents a number that provides all common
 *              mathmatical operations.
 *
 * @since       0.1.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        30.05.15
 */
template< typename T >
class
vector< T, if_pod_type< T > >
{
public:
    // type declarations
    typedef T                                          pod_type; //!< The POD type of matrix elements
    typedef typename smart_array< pod_type >::iterator iterator; //!< Vector iterator
    
    /*!
     * @brief       Specifing the type a vector can have
     * @details     Indicator for the type of the current vector.
     *              A vector can either be a column vector or a
     *              row vector.
     *
     * @ingroup     vector
     */
    enum type
    {
        ROW,        //!< Represents a column vector of dimension \f$M\times 1\f$ where \f$M\in\mathbb{N}^+\f$
        COLUMN      //!< Represents a row vector of dimension \f$1\times M\f$ where \f$M\in\mathbb{N}^+\f$
    };
    
    // ivars
    const smart_array< pod_type > mem;  //!< vector data
    const size_t size;                  //!< size of vector
    const type type;                    //!< type of vector
    
    // methods
    inline                                      vector();
    inline                                      vector(const size_t& s, const enum type& t = vector< pod_type >::ROW);
    inline                                      vector(const vector< pod_type >& vec);
    inline                                      vector(const vector< pod_type >& vec, const enum type& t);
    inline                                     ~vector();
    
    inline       matrix< pod_type >             operator*(const vector< pod_type >& v);
    inline       vector< pod_type >             operator*(const pod_type& s);
    inline       vector< complex< pod_type > >  operator*(const complex< pod_type >& s);
    
    inline       vector< pod_type >             operator+();
    inline       vector< pod_type >             operator-();
    
    inline       vector< complex< pod_type > >  operator*(const matrix< complex< pod_type > >& mat);
    
    inline const vector< pod_type >&            operator=(const vector< pod_type >& v);
    
    inline const vector< pod_type >&            operator*=(const vector< pod_type >& v);
    inline const vector< pod_type >&            operator*=(const pod_type& s);
    
    inline       pod_type&                      operator[](const size_t& idx);
    inline const pod_type&                      operator[](const size_t& idx) const;
    
    inline       pod_type&                      operator()(const size_t& idx);
    inline const pod_type&                      operator()(const size_t& idx) const;
    
    inline       void                           transpose();
    
    inline       iterator                       begin();
    inline       iterator                       end();
};

/*!
 * @brief           Default constructor for a vector
 * @details         Constructs an empty vector which represents a row vector.
 *                  The memory is not initialized and can contain random
 *                  contents!
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::vector()
    : type(vector< pod_type >::ROW)
    , size(0)
{}

/*!
 * @brief           Constructor for creating a vector of size \f$s\f$
 * @details         Constructs a vector of type \f$t\f$ and size \f$s\f$.
 *                  The memory is not initialized and can contain random
 *                  contents!
 *
 * @param[in]       s The size of the created vector
 * @param[in]       type Type of the created vector where \f$t\f$ can be
 *                  either type::ROW or type::COLUMN
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::vector(const size_t& s, const enum type& type)
    : size(s)
    , type(type)
{
    access::rw(mem) = smart_array< pod_type >(s);
}

/*!
 * @brief           The copy constructor for the vector class.
 * @details         Constructs a new vector that is a copy of the given
 *                  vector. The newly created vector is in the same state
 *                  than the given vector.
 *
 * @param[in]       vec The vector that is supposed to be copied.
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::vector(const vector< pod_type >& vec)
    : size(vec.size)
    , type(vec.type)
{
    access::rw(mem) = vec.mem;
}

/*!
 * @brief           The copy constructor for a vector and a new type.
 * @details         Constructs a new vector that is a copy of the given
 *                  vector but with a given new type \f$t\f$. The newly
 *                  created vector is in the same state than the given
 *                  vector.
 *
 * @param[in]       vec The vector that is supposed to be copied.
 * @param[in]       type The new type of the new vector which can either be
 *                  type::ROW or type::COlUMN.
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::vector(const vector< pod_type >& vec, const enum type& type)
    : size(vec.size)
    , type(type)
{
    access::rw(mem) = vec.mem;
}

/*!
 * @brief           The destructor for the vector.
 * @details         Frees the allocated memory of the vector and preparing
 *                  the object to be released safely.
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::~vector()
{}

/*!
 * @brief           Multiplication operator for a vector and a scalar.
 * @details         Multplies a given scalar to each element of the current
 *                  vector.
 *
 * @param[in]       s The given scalar value.
 *
 * @return          A new vector containing the product of each vector element
 *                  and the scalar.
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator*(const pod_type& s)
{
    vector< pod_type > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] * s;
    }
    
    return result;
}

/*!
 * @brief           Multiplication operator for a vector and a complex scalar.
 * @details         Multiplies the given complex scalar to each element
 *                  in the current vector.
 *
 * @param[in]       s The given complex scalar.
 *
 * @return          A complex vector containing the product of each
 *                  vector element and the scalar value.
 */
template< typename T >
inline
vector< complex< T > > vector< T, if_pod_type< T > >::operator*(const complex< pod_type >& s)
{
    vector< complex< pod_type > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = complex< pod_type >(mem[i], 0) * s;
    }
    
    return result;
}



/*!
 * @brief           Positive sign operator for a vector.
 * @details         returns the vector as it is.
 *
 * @return          returns the vector as it is
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator+()
{
    vector< pod_type > result;
    access::rw(result.mem) = mem;
    
    return result;
}

/*!
 * @brief           Negative sign operator for a vector.
 * @details         Negates each element of the current vector.
 *
 * @return          A vector where each element is negated.
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator-()
{
    vector< pod_type > result;
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = -mem[i];
    }
    
    return result;
}



/*!
 * @brief           The assignment operator for a vector.
 * @details         Assigns the given vector to the current vector.
 *                  The newly created vector is in the same state than
 *                  the given vector.
 *
 * @param[in]       v The given vector that is supposed to be assignend.
 *
 * @return          A reference to the current vector.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator=(const vector< pod_type >& v)
{
    if ( this == &v )
    {
        return *this;
    }
    
    access::rw(size) = v.size;
    access::rw(type) = v.type;
    access::rw(mem)  = v.mem;
    
    return *this;
}

/*!
 * @brief           The multiplication assignment operator for a vector and a scalar.
 * @details         Multiplies a given scalar to the current vector. The scalar value
 *                  gets multiplied to each vector element.
 *
 * @param[in]       s The scalar that is supposed to be multiplied to the current vector.
 *
 * @return          The reference to the current vector.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator*=(const pod_type& s)
{
    size_t i;
    for (i = 0; i < size; ++i)
    {
        access::rw(mem[i]) *= s;
    }
    
    return *this;
}

/*!
 * @brief           The write subscript operator for the vector object.
 * @details         Overload for the subscript operator. The subscript operator is used
 *                  to access the elements of the vector by using square brackets on the
 *                  vector object like this.
 *                  @code
 *                      vector<double> v(5, vector< double >::ROW);
 *                      v.fill(5);
 *
 *                      // use subscript to overwrite the third value
 *                      v[3] = 7;
 *
 *                      // The content of v looks like this
 *                      // v =
 *                      //     5  5  7  5  5
 *                  @endcode
 *
 * @param[in]       idx The index of the element that is supposed to be returned.
 *
 * @return          The reference to element at index idx
 */
template< typename T >
inline
T& vector< T, if_pod_type< T > >::operator[](const size_t& idx)
{
    return access::rw(mem[idx]);
}

/*!
 * @brief           The read subscript operator for the vector objects.
 * @details         Overload for the subscript operator. The subscript operator is used
 *                  to access the elements of the vector by using square brackets on the
 *                  vector object like this.
 *                  @code
 *                      double a = 3;
 *                      vector<double> v(5, vector<double>::ROW)
 *                      v.fill(5);
 *
 *                      // use subscript to read from the vector
 *                      a = v[5];
 *
 *                      // The content of a looks like this
 *                      // a = 5
 *                  @endcode
 *
 * @param[in]       idx The index of the element that is supposed to be read.
 *
 * @return          The element at index idx
 */
template< typename T >
inline
const T& vector< T, if_pod_type< T > >::operator[](const size_t& idx) const
{
    return access::rw(mem[idx]);
}



/*!
 * @brief           The write access operator for the vector object.
 * @details         Overload for the access operator. The access operator is used
 *                  to access the elements of the vector by using parenthesis on the
 *                  vector object like this.
 *                  @code
 *                      vector<double> v(5, vector<double>::ROW);
 *                      v.fill(5);
 *
 *                      // use subscript to overwrite the third value
 *                      v(3) = 7;
 *
 *                      // The content of v looks like this
 *                      // v =
 *                      //     5  5  7  5  5
 *                  @endcode
 *
 * @param[in]       idx The index of the element that is supposed to be returned.
 *
 * @return          The reference to element at index idx
 */
template< typename T >
inline
T& vector< T, if_pod_type< T > >::operator()(const size_t& idx)
{
    return mem[idx];
}

/*!
 * @brief           The read access operator for the vector objects.
 * @details         Overload for the access operator. The access operator is used
 *                  to access the elements of the vector by using parenthesis on the
 *                  vector object like this.
 *                  @code
 *                      double a = 3;
 *                      vector<double> v(5, vector<double>::ROW)
 *                      v.fill(5);
 *
 *                      // use subscript to read from the vector
 *                      a = v(5);
 *
 *                      // The content of a looks like this
 *                      // a = 5
 *                  @endcode
 *
 * @param[in]       idx The index of the element that is supposed to be read.
 *
 * @return          The element at index idx
 */
template< typename T >
inline
const T& vector< T, if_pod_type< T > >::operator()(const size_t& idx) const
{
    return mem[idx];
}

/*!
 * @brief           Transposes the current vector.
 * @details         Switching vector type from type::ROW to type::COLUMN or
 *                  from type::COLUMN to type::ROW
 */
template< typename T >
inline
void vector< T, if_pod_type< T > >::transpose()
{
    type = (type == vector< pod_type >::ROW) ? vector< pod_type >::COLUMN : vector< pod_type >::ROW;
}

template< typename T >
inline
typename vector< T, if_pod_type< T > >::iterator vector< T, if_pod_type< T > >::begin()
{
    return access::rw(mem).begin();
}

template< typename T >
inline
typename vector< T, if_pod_type< T > >::iterator vector< T, if_pod_type< T > >::end()
{
    return access::rw(mem).end();
}


/*!
 * @brief           Outstream operator overload.
 * @details         The out-steam operator is used to print the vector
 *                  in a nice form over the std::cout stream.
 *
 * @param[in,out]   o The stream object
 * @param[in]       v The vector that should be printed
 *
 * @return          The reference to the given out-stream.
 */
template< typename S >
std::ostream& operator<<(std::ostream& o, const vector< S >& v)
{
    // setting decimal precesion
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl;
    
    int width   = 10;
    auto format = std::fixed;
    
    // reduce size for integers
    if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
    {
        width = 5;
    }
    
    // check values
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        S val = v[i];
        if (std::abs(val) >= 10)
        {
            width   = 11;
            format  = std::fixed;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 6;
            }
        }
        
        if (std::abs(val) >= 100)
        {
            width   = 12;
            format  = std::fixed;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 7;
            }
        }
        
        if (std::abs(val) >= 1000)
        {
            width   = 14;
            format  = std::scientific;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 10;
            }
        }
    }
    
    if (v.type == vector< S >::ROW)
    {
        for (i = 0; i < v.size; ++i)
        {
            // get entry
            S val = v[i];
            
            // create string
            o << std::setfill(' ');
            o << std::right << std::setw(width);
            o << format << std::setprecision(4) << val;
        }
        o << std::endl;
    }
    else
    {
        for (i = 0; i < v.size; ++i)
        {
            // get entry
            S val = v[i];
            
            // create string
            o << std::setfill(' ');
            o << std::right << std::setw(width);
            o << format << std::setprecision(4) << val << std::endl;
        }
    }
    
    std::cout.flags( f );
    return o;
}

/*!
 * @}
 */

PSOFFT_END

#endif /* vector.hpp */
