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

#ifndef PSOFFTlib_vector_cx_hpp
#define PSOFFTlib_vector_cx_hpp

PSOFFT_BEGIN

/*!
 * @ingroup    vector
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
vector< complex< T >, if_pod_type< T > >
{
public:
    // type declarations
    typedef T pod_type;                 //!< The POD type of matrix elements
    
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
    
    // iterators
    typedef typename smart_array< complex< pod_type > >::iterator iterator; //!< Vector iterator
    
    // ivars
    const smart_array< complex< pod_type > > mem;   //!< vector data
    const size_t size;                              //!< size of vector
    const type type;                                //!< type of vector
    
    inline                                vector();
    inline                                vector(const size_t& s, const enum type& type = vector< complex< T > >::ROW);
    inline                                vector(const vector< T >& vec);
    inline                                vector(const vector< complex< T > >& vec);
    inline                                vector(const vector< T >& vec, const enum type& type);
    inline                                vector(const vector< complex< T > >& vec, const enum type& type);
    inline                                vector(const size_t& s, const pod_type& initial, const enum type& type = vector< complex< pod_type > >::ROW);
    inline                               ~vector();
    
    inline       vector< complex< T > >   operator*(const T& s);
    inline       vector< complex< T > >   operator*(const complex< T >& s);
    
    inline       vector< complex< T > >   operator+();
    inline       vector< complex< T > >   operator-();
    
    inline       vector< complex< T > >   operator*(const matrix< T >& mat);
    inline       vector< complex< T > >   operator*(const matrix< complex< T > >& mat);
    
    inline const vector< complex< T > >&  operator=(const vector< T >& v);
    inline const vector< complex< T > >&  operator=(const vector< complex< T > >& v);
    
    inline const vector< complex< T > >&  operator*=(const T& s);
    
    inline       vector< complex< T > >   operator-(const vector< complex< T > >& v);
    
    inline       complex< T >&            operator[](const size_t& idx);
    inline const complex< T >&            operator[](const size_t& idx) const;
    
    inline       complex< T >&            operator()(const size_t& idx);
    inline const complex< T >&            operator()(const size_t& idx) const;
    
    inline       void                     transpose();
    
    // iterator methods
    inline       iterator                 begin();
    inline       iterator                 end();
};

template< typename T >
inline
vector< complex< T >, if_pod_type< T > >::vector()
    : size(0)
    , type(vector< complex< pod_type > >::ROW)
{}

template< typename T >
inline
vector< complex< T >, if_pod_type< T > >::vector(const size_t& s, const enum type& type)
    : size(s)
    , type(type)
{
    access::rw(mem) = smart_array< complex< pod_type > >(s);
}

template< typename T >
inline
vector< complex< T >, if_pod_type< T > >::vector(const vector< T >& vec)
    : size(vec.size)
    , type(vec.type)
{
    access::rw(mem) = smart_array< complex< pod_type > >(vec.size);
    
    size_t i;
    for (i = 0; i < vec.size; ++i)
    {
        access::rw(mem[i]) = complex< T >(vec[i], 0);
    }
}

template< typename T >
inline
vector< complex< T >, if_pod_type< T > >::vector(const vector< complex< T > >& vec)
    : size(vec.size)
    , type(vec.type)
{
    access::rw(mem) = vec.mem;
}

template< typename T >
inline
vector< complex< T >, if_pod_type< T > >::vector(const vector< T >& vec, const enum type& type)
    : size(vec.size)
    , type(type)
{
    access::rw(mem) = smart_array< complex< pod_type > >(vec.size);
    
    size_t i;
    for (i = 0; i < vec.size; ++i)
    {
        access::rw(mem[i]) = complex< T >(vec[i], 0);
    }
}

template< typename T >
inline
vector< complex< T >, if_pod_type< T > >::vector(const vector< complex< T > >& vec, const enum type& type)
    : size(vec.size)
    , type(type)
{
    access::rw(mem) = vec.mem;
}

template< typename T >
inline
vector< complex< T >, if_pod_type< T > >::vector(const size_t& s, const T& initial, const enum type& type)
    : size(s)
    , type(type)
{
    access::rw(mem) = smart_array< complex< pod_type > >(size);
    
    if (size > 0)
    {
        for (typename smart_array< complex< pod_type > >::iterator it = access::rw(mem).begin(); it != access::rw(mem).end(); ++it)
        {
            access::rw(*it) = 0;
        }
    }
}

template< typename T >
inline
vector< complex< T >, if_pod_type< T > >::~vector()
{}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, if_pod_type< T > >::operator*(const T& s)
{
    vector< complex< pod_type > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = access::rw(mem[i]) * complex< pod_type >(s, 0);
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, if_pod_type< T > >::operator*(const complex< T >& s)
{
    vector< complex< pod_type > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] * s;
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, if_pod_type< T > >::operator+()
{
    vector< complex< pod_type > > result(size, type);
    access::rw(result.mem) = mem;
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, if_pod_type< T > >::operator-()
{
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = - mem[i];
    }
    
    return result;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, if_pod_type< T > >::operator=(const vector< T >& v)
{
    size = v.size;
    type = v.type;
    
    mem = smart_array< complex< pod_type > >(size);
    
    if (size > 0)
    {
        size_t i;
        for (i = 0; i < size; ++i)
        {
            mem[i] = complex< T >(v[i], 0);
        }
    }
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, if_pod_type< T > >::operator=(const vector< complex< T > >& v)
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

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, if_pod_type< T > >::operator*=(const T& s)
{
    
    for (typename smart_array< complex< pod_type > >::iterator it = access::rw(mem).begin(); it != access::rw(mem).end(); ++it)
    {
        *it *= complex< pod_type >(s, 0);
    }
    
    return *this;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, if_pod_type< T > >::operator-(const vector< complex< T > >& v)
{
    if ( size != v.size || type != v.type)
    {
        psofft_error("%s", "size mismatch in complex vector-vector subtraction.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = access::rw(v[i]) - access::rw(mem[i]);
    }
    
    return result;
}

template< typename T >
inline
complex< T >& vector< complex< T >, if_pod_type< T > >::operator[](const size_t& idx)
{
    return access::rw(mem[idx]);
}

template< typename T >
inline
const complex< T >& vector< complex< T >, if_pod_type< T > >::operator[](const size_t& idx) const
{
    return access::rw(mem[idx]);
}

template< typename T >
inline
complex< T >& vector< complex< T >, if_pod_type< T > >::operator()(const size_t& idx)
{
    return mem[idx];
}

template< typename T >
inline
void vector< complex< T >, if_pod_type< T > >::transpose()
{
    type = (type == vector< complex< T > >::ROW) ? vector< complex< T > >::COLUMN : vector< complex< T > >::ROW;
}

template< typename T >
inline
typename vector< complex< T >, if_pod_type< T > >::iterator vector< complex< T >, if_pod_type< T > >::begin()
{
    return access::rw(mem).begin();
}

template< typename T >
inline
typename vector< complex< T >, if_pod_type< T > >::iterator vector< complex< T >, if_pod_type< T > >::end()
{
    return access::rw(mem).end();
}


template< typename S >
std::ostream& operator<<(std::ostream& o, const vector< complex< S > >& v)
{
    // setting decimal precesion
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl;
    
    int width   = 20;
    auto format = std::fixed;
    
    if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
    {
        width = 10;
    }
    
    // check values
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        complex< S > c = v[i];
        
        if (std::abs(c.re) >= 10 || std::abs(c.im) >= 10)
        {
            width   = 22;
            format  = std::fixed;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 12;
            }
        }
        
        if (std::abs(c.re) >= 100 || std::abs(c.im) >= 100)
        {
            width   = 24;
            format  = std::fixed;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 14;
            }
        }
        
        if (std::abs(c.re) >= 1000 || std::abs(c.im) >= 1000)
        {
            width   = 28;
            format  = std::scientific;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 18;
            }
        }
    }
    
    // prepare output and print
    if (v.type == vector< complex< S > >::ROW)
    {
        for (i = 0; i < v.size; ++i)
        {
            // get entry
            complex< S > c = v[i];
            
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
        o << std::endl;
    }
    else
    {
        for (i = 0; i < v.size; ++i)
        {
            // get entry
            complex< S > c = v[i];
            
            // create string
            std::ostringstream val;
            
            // add real value to string
            val << format << std::setprecision(4) << c.re;
            val << (c.im < 0 ? " - " : " + ") << (c.im == 0 ?  0 : std::abs(c.im)) << "i";
            
            // get string from stream
            std::string str = val.str();
            
            // set filling character
            o << std::setfill(' ') << std::right << std::setw(width) << str << std::endl;
        }
    }
    
    std::cout.flags( f );
    return o;
}

/*!
 * @}
 */

PSOFFT_END

#endif /* cx_vector_dec.hpp */
