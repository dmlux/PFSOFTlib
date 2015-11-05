//
//  random.hpp
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

#ifndef PDSOFTlib_random_hpp
#define PDSOFTlib_random_hpp

PDSOFT_BEGIN

/*!
 * @brief       A context for random number creation
 */
template< typename pod_type, typename derived >
struct randctx
{
    inline const derived& get_ref() const;
};

template< typename pod_type, typename derived >
inline
const derived& randctx< pod_type, derived >::get_ref() const
{
    return static_cast< const derived& >(*this);
}

template< typename T >
struct uniform_real_distribution< T, typename if_true< is_real_type< T >::value >::type > : randctx< T, uniform_real_distribution< T > >
{
public:
    random_engine engine = random_engine::DEFAULT;  //!< Random engine that should be used
    T min;                                          //!< Minimum random value
    T max;                                          //!< Maximum random value
};

template< typename T >
struct uniform_int_distribution< T, typename if_true< is_integral_type< T >::value >::type > : randctx< T, uniform_int_distribution< T > >
{
public:
    random_engine engine = random_engine::DEFAULT;  //!< Random engine that should be used
    T min;                                          //!< Minimum random value
    T max;                                          //!< Maximum random value
};


/*!
 * @brief           filling Fourier coefficients container with random values
 * @details         Fills the fourier coefficients container with random values
 *                  in the given range.
 *
 * @param[in, out]  fc The fourier coefficients container
 * @param[in]       min The minimum value of the random number range
 * @param[in]       max The maximum value of the random number range
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            18.06.15
 *
 * @since           0.1.1
 */
template< typename pod_type, typename distribution >
inline
uniform_real_dist_type< distribution > rand(DSOFTFourierCoefficients& fc, const randctx< pod_type, distribution >& ctx)
{
    // cast to correct underlying type
    uniform_real_distribution< pod_type > dist = ctx.get_ref();
    
    // get min and max values
    pod_type min = dist.min < dist.max ? dist.min : dist.max;
    pod_type max = dist.min < dist.max ? dist.max : dist.min;
    
    // create timeval objects
    struct timeval tv;
    
    // get current time in microseconds
    gettimeofday(&tv, NULL);
    
    // create seed
    unsigned long seed = 1000000L * tv.tv_sec + tv.tv_usec;
    
    // C++11 random numbers uniformly distributed
    if ( dist.engine == random_engine::DEFAULT || dist.engine == random_engine::MERSENNE_TWISTER64 )
    {
        std::mt19937_64 e(seed);
        std::uniform_real_distribution< pod_type > d(min, max);
        
        // iterate over degree
        for (int l = 0; l < fc.bandwidth; ++l)
        {
            // iterate over M order
            for (int M = -l; M <= l; ++M)
            {
                // iterate over M' order
                for (int Mp = -l; Mp <= l; ++Mp)
                {
                    fc(l,M,Mp).re = d(e);
                    fc(l,M,Mp).im = d(e);
                }
            }
        }
    }
    
    else if (dist.engine == random_engine::MERSENNE_TWISTER)
    {
        std::mt19937 e(seed);
        std::uniform_real_distribution< pod_type > d(min, max);
        
        // iterate over degree
        for (int l = 0; l < fc.bandwidth; ++l)
        {
            // iterate over M order
            for (int M = -l; M <= l; ++M)
            {
                // iterate over M' order
                for (int Mp = -l; Mp <= l; ++Mp)
                {
                    fc(l,M,Mp).re = d(e);
                    fc(l,M,Mp).im = d(e);
                }
            }
        }
    }
}

/*!
 * @brief           Fills a given complex vector with real random values.
 * @details         The given vector gets filled with random complex values in
 *                  range of \f$[\min, \max]\f$
 *
 * @param[in, out]  vec The complex vector that is supposed to be filled with random
 *                  values.
 * @param[in]       min The min value for the random co-domain.
 * @param[in]       max The max value for the random co-domain.
 * @tparam          T The element type of the vector. The type has to be a floating
 *                  point type (float, double or long double).
 *
 * @note            Only available for vectors of type float, double or long double!
 *
 * @ingroup         vector
 */
template< typename T >
inline
void_real_type< T > rand(vector< complex< T > >& vec, const double& min, const double& max)
{
    if (min > max)
    {
        pdsoft_warning("%s", "min value is greater than max value in rand function for complex vectors.");
        return;
    }
    
    // create timeval object
    struct timeval tv;
    
    // get current time in microseconds
    gettimeofday(&tv, NULL);
    
    // create seed
    unsigned long seed = 1000000L * tv.tv_sec + tv.tv_usec;
    
    // C++11 random numbers uniformly distributed
    std::default_random_engine e(seed);
    std::uniform_real_distribution< T > d(min, max);
    
    // fill vector with randoms
    size_t i;
    for (i = 0; i < vec.size; ++i)
    {
        vec[i].re = d(e);
        vec[i].im = d(e);
    }
}

PDSOFT_END

#endif /* random.hpp */
