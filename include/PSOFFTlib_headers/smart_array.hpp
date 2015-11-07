//
//  smart_array.hpp
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

#ifndef PSOFFTlib_smart_array_hpp
#define PSOFFTlib_smart_array_hpp

PSOFFT_BEGIN

template< typename T >
class smart_array
{
public:
    
    // typedefs
    typedef T         pod_type;
    typedef pod_type* iterator;
    
private:
    // properties of the array
    pod_type* memory;   //!< internal storage
    size_t size;        //!< size of this array
    
public:
    // methods
    inline                                 smart_array();
    inline                                 smart_array(const size_t& size);
    inline                                 smart_array(const smart_array< pod_type >& array);
    inline                                 smart_array(smart_array< pod_type >&& array);
    inline                                ~smart_array();
    
    inline       pod_type&                 operator*();
    inline       pod_type*&                operator+(const size_t& idx);
    inline       pod_type*&                operator-(const size_t& idx);
    
    inline       pod_type&                 operator[](const size_t& idx);
    inline const pod_type&                 operator[](const size_t& idx) const;
    
    inline const smart_array< pod_type >&  operator=(const smart_array< pod_type >& rhs);
    inline const smart_array< pod_type >&  operator=(smart_array< pod_type >&& rhs);
    
    inline       iterator                  begin();
    inline       iterator                  end();
};

/*!
 * @brief           Constructor for a smart array of no size
 */
template< typename T >
inline
smart_array< T >::smart_array()
    : memory(nullptr)
    , size(0)
{}

/*!
 * @brief           Constructor for a smart array of given size
 */
template< typename T >
inline
smart_array< T >::smart_array(const size_t& s)
{
    memory = (pod_type*) malloc(s * sizeof(pod_type));
    size   = s;
}

/*
 * @brief           Constructor for a smart array with same contents
 *                  and properties as another given array.
 */
template< typename T >
inline
smart_array< T >::smart_array(const smart_array< pod_type >& array)
{
    memory = (pod_type*) malloc(array.size * sizeof(pod_type));
    size   = array.size;
    
    memcpy(memory, array.memory, size * sizeof(pod_type));
}

/*!
 * 
 */
template< typename T >
inline
smart_array< T >::smart_array(smart_array< pod_type >&& array)
{
    pod_type* tmp = memory;
    memory        = array.memory;
    array.memory  = tmp;
}

template< typename T >
inline
smart_array< T >::~smart_array()
{
    free( memory );
}

template< typename T >
inline
T& smart_array< T >::operator*()
{
    return *memory;
}

template< typename T >
inline
T*& smart_array< T >::operator+(const size_t& idx)
{
    return memory + idx;
}

template< typename T >
inline
T*& smart_array< T >::operator-(const size_t& idx)
{
    return memory - idx;
}

template< typename T >
inline
T& smart_array< T >::operator[](const size_t& idx)
{
    return *(memory + idx);
}

template< typename T >
inline
const T& smart_array< T >::operator[](const size_t& idx) const
{
    return *(memory + idx);
}

template< typename T >
inline
const smart_array< T >& smart_array< T >::operator=(const smart_array< pod_type >& rhs)
{
    free( memory );
    
    size   = rhs.size;
    memory = (pod_type*) malloc(size * sizeof(pod_type));
    
    memcpy(memory, rhs.memory, size * sizeof(pod_type));
    
    return *this;
}

template< typename T >
inline
const smart_array< T >& smart_array< T >::operator=(smart_array< pod_type >&& rhs)
{
    pod_type* tmp = memory;
    memory        = rhs.memory;
    rhs.memory    = tmp;
    size          = rhs.size;
    
    return *this;
}


template< typename T >
inline
typename smart_array< T >::iterator smart_array< T >::begin()
{
    return memory;
}

template< typename T >
inline
typename smart_array< T >::iterator smart_array< T >::end()
{
    return memory + size;
}

PSOFFT_END

#endif /* smart_array.hpp */
