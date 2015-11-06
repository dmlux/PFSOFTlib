//
//  forward_declarations.hpp
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

#ifndef PSOFFTlib_forward_declarations_hpp
#define PSOFFTlib_forward_declarations_hpp

PSOFFT_BEGIN

// check if type is num type
template< typename T > struct is_num_type                               { static const bool value = false;  };
template<>             struct is_num_type< short >                      { static const bool value = true;   };
template<>             struct is_num_type< int >                        { static const bool value = true;   };
template<>             struct is_num_type< long >                       { static const bool value = true;   };
template<>             struct is_num_type< long long >                  { static const bool value = true;   };
template<>             struct is_num_type< unsigned short >             { static const bool value = true;   };
template<>             struct is_num_type< unsigned int >               { static const bool value = true;   };
template<>             struct is_num_type< unsigned long >              { static const bool value = true;   };
template<>             struct is_num_type< unsigned long long >         { static const bool value = true;   };
template<>             struct is_num_type< float >                      { static const bool value = true;   };
template<>             struct is_num_type< double >                     { static const bool value = true;   };
template<>             struct is_num_type< long double >                { static const bool value = true;   };

// check if template is true
template< bool T, typename U = void > struct if_true                    {                                   };
template< typename U >                struct if_true< true, U >         { typedef U type;                   };

// template alias
template< typename pod > using if_pod_type = typename if_true< is_num_type< pod >::value >::type;

// forward declarations
template< typename >                        class  complex;

// matrix
template< typename, typename = void >       class  matrix;
template< typename T >                      class  matrix< T,            if_pod_type< T > >;
template< typename T >                      class  matrix< complex< T >, if_pod_type< T > >;

// vector
template< typename, typename = void >       class  vector;
template< typename T >                      class  vector< T,            if_pod_type< T > >;
template< typename T >                      class  vector< complex< T >, if_pod_type< T > >;

// 3D grid
template< typename, typename = void >       struct grid3D;
template< typename T >                      struct grid3D< complex< T >, if_pod_type< T > >;

// smart_array
template< typename T >                      class  smart_array;
                                            class  stopwatch;

                                            struct DSOFTFourierCoefficients;

template< typename, typename >              struct randctx;
template< typename, typename = void >       struct uniform_int_distribution;
template< typename, typename = void >       struct uniform_real_distribution;
template< typename, typename = void >       struct normal_distribution;

/*!
 * @brief       A collection of usable random engines
 */
enum class random_engine
{
    DEFAULT,            //!< Mersenne Twister 64 bit generator
    MINSTD_RAND,        //!< Minimal Standard generator
    MINSTD_RAND0,       //!< Minimal Standard 0 generator
    MERSENNE_TWISTER,   //!< Mersenne Twister generator
    MERSENNE_TWISTER64, //!< Mersenne Twister 64 bit generator
    RANLUX24_BASE,      //!< Ranlux 24 Base generator
    RANLUX48_BASE,      //!< Ranlux 48 Base generator
    RANLUX24,           //!< Ranlux 24 generator
    RANLUX48,           //!< Ranlux 48 generator
    KNUTH_B             //!< Knuth B generator
};

/*!
 * @brief       A constant variable accessor
 */
class access
{
public:
    
    /*!
     * @brief           Removes the constantness of given element
     */
    template< typename T > psofft_inline static T&  rw (const T& x)        { return const_cast< T&  >(x); }
    template< typename T > psofft_inline static T*& rwp(const T* const& x) { return const_cast< T*& >(x); }
};

/*!
 * @brief       Collection of trait classes.
 * @details     Trait classes can be used to detect the type of template
 *              variables at compile time to optimize function efficiency.
 * @defgroup    traits Trait structs
 * @{
 */

// check for same type
template< typename T, typename U > struct same_type                         { static const bool value = false;  };
template< typename T >             struct same_type< T, T >                 { static const bool value = true;   };

// check for different types
template< typename T, typename U > struct different_type                    { static const bool value = true;   };
template< typename T >             struct different_type< T, T >            { static const bool value = false;  };

// check if type is complex number
template< typename T > struct is_complex                                    { static const bool value = false;  };
template<>             struct is_complex< complex< short > >                { static const bool value = true;   };
template<>             struct is_complex< complex< int > >                  { static const bool value = true;   };
template<>             struct is_complex< complex< long > >                 { static const bool value = true;   };
template<>             struct is_complex< complex< long long > >            { static const bool value = true;   };
template<>             struct is_complex< complex< unsigned short > >       { static const bool value = true;   };
template<>             struct is_complex< complex< unsigned int > >         { static const bool value = true;   };
template<>             struct is_complex< complex< unsigned long > >        { static const bool value = true;   };
template<>             struct is_complex< complex< unsigned long long > >   { static const bool value = true;   };
template<>             struct is_complex< complex< float > >                { static const bool value = true;   };
template<>             struct is_complex< complex< double > >               { static const bool value = true;   };
template<>             struct is_complex< complex< long double > >          { static const bool value = true;   };

// check if type is real
template< typename T > struct is_real_type                                  { static const bool value = false;  };
template<>             struct is_real_type< float >                         { static const bool value = true;   };
template<>             struct is_real_type< double >                        { static const bool value = true;   };
template<>             struct is_real_type< long double >                   { static const bool value = true;   };

// check if type is integral
template< typename T > struct is_integral_type                              { static const bool value = false;  };
template<>             struct is_integral_type< short >                     { static const bool value = true;   };
template<>             struct is_integral_type< int >                       { static const bool value = true;   };
template<>             struct is_integral_type< long >                      { static const bool value = true;   };
template<>             struct is_integral_type< long long >                 { static const bool value = true;   };
template<>             struct is_integral_type< unsigned short >            { static const bool value = true;   };
template<>             struct is_integral_type< unsigned int >              { static const bool value = true;   };
template<>             struct is_integral_type< unsigned long >             { static const bool value = true;   };
template<>             struct is_integral_type< unsigned long long >        { static const bool value = true;   };

/*!
 * @}
 */

/*!
 * @brief       Collection of restrictors to restrict element types for special functions
 * @defgroup    restrictors Restrictors
 * @{
 */

// Real values for return type void
template< typename T > struct void_real_only                                    {                                               };
template<>             struct void_real_only< float >                           { typedef void result;                          };
template<>             struct void_real_only< double >                          { typedef void result;                          };
template<>             struct void_real_only< long double >                     { typedef void result;                          };

template< typename T > using void_real_type = typename void_real_only< T >::result;

// Numbers only
template< typename T > struct numbers_only                                      {                                               };
template<>             struct numbers_only< short >                             { typedef short result;                         };
template<>             struct numbers_only< int >                               { typedef int result;                           };
template<>             struct numbers_only< long >                              { typedef long result;                          };
template<>             struct numbers_only< long long >                         { typedef long long result;                     };
template<>             struct numbers_only< unsigned short >                    { typedef unsigned short result;                };
template<>             struct numbers_only< unsigned int >                      { typedef unsigned int result;                  };
template<>             struct numbers_only< unsigned long >                     { typedef unsigned long result;                 };
template<>             struct numbers_only< unsigned long long >                { typedef unsigned long long result;            };
template<>             struct numbers_only< float >                             { typedef float result;                         };
template<>             struct numbers_only< double >                            { typedef double result;                        };
template<>             struct numbers_only< long double >                       { typedef long double result;                   };

template< typename T > using number_type = typename numbers_only< T >::result;

// Real numbers for return type void
template< typename T > struct void_numbers_only                                 {                                               };
template<>             struct void_numbers_only< short >                        { typedef void result;                          };
template<>             struct void_numbers_only< int >                          { typedef void result;                          };
template<>             struct void_numbers_only< long >                         { typedef void result;                          };
template<>             struct void_numbers_only< long long >                    { typedef void result;                          };
template<>             struct void_numbers_only< unsigned short >               { typedef void result;                          };
template<>             struct void_numbers_only< unsigned int >                 { typedef void result;                          };
template<>             struct void_numbers_only< unsigned long >                { typedef void result;                          };
template<>             struct void_numbers_only< unsigned long long >           { typedef void result;                          };
template<>             struct void_numbers_only< float >                        { typedef void result;                          };
template<>             struct void_numbers_only< double >                       { typedef void result;                          };
template<>             struct void_numbers_only< long double >                  { typedef void result;                          };

template< typename T > using void_number_type = typename void_numbers_only< T >::result;

// Uniform real dist types
template< typename T > struct uniform_real_dist_only                                                {                           };
template<>             struct uniform_real_dist_only< uniform_real_distribution< float > >          { typedef void result;      };
template<>             struct uniform_real_dist_only< uniform_real_distribution< double > >         { typedef void result;      };
template<>             struct uniform_real_dist_only< uniform_real_distribution< long double > >    { typedef void result;      };

template< typename T > using uniform_real_dist_type = typename uniform_real_dist_only< T >::result;

/*!
 * @}
 */

/*!
 * @brief       A collection of usable constants
 */
template< typename T >
class constants
{
public:
    static const T pi; //!< The ratio of a circle's circumference to its diameter
    static const T e;  //!< e is the limit of (1 + 1/n)^n for n to infinity
};

template< typename T >
const T constants< T >::pi = static_cast< T >(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664);

template< typename T >
const T constants< T >::e = static_cast< T >(2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992);

PSOFFT_END

#endif /* forward_declarations.hpp */
