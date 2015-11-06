//
//  comlex.hpp
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

#ifndef PSOFFTlib_complex_hpp
#define PSOFFTlib_complex_hpp

PSOFFT_BEGIN

/*!
 * @brief       Collection of classes and functions for complex numbers for
 *              mathematical purposes.
 * @details     Contains a set of classes and functions that are useful for 
 *              dealing with complex numbers.
 * @defgroup    complex Complex
 * @{
 */

/*! 
 * @brief   A C++ implementation of a complex number
 * @details A complex number that contains all basic operations for calculating 
 *          with this type of numbers. All necessary operators are overloaded to 
 *          provide easier readability in algorihtms. The complex number uses a 
 *          template parameter to store the real and imaginary part in a member 
 *          variable of a specific type. Any POD type of C or C++ can be used to 
 *          store those parts of the complex number. Additionally to the build in 
 *          operators the calculation between PODs and a complex number are 
 *          implemented to.
 *
 * @tparam  T An element type which represents a number that provides all common
 *          mathmatical operations.
 *
 * @since   0.0.1
 *
 * @todo    Add additional operators for different datatypes.
 * @todo    Write conversion operator overloads.
 *
 * @author  Denis-Michael Lux <denis.lux@icloud.com>
 * @date    11.01.15
 */
template< typename T >
class complex
{
public:
    typedef T pod_type;
    
    pod_type re;   //!< The real part of the complex number
    pod_type im;   //!< The imaginary part of the complex number
    
    inline                            ~complex();
    inline                             complex();
    inline                             complex(const pod_type real_imag);
    inline                             complex(const pod_type real, const pod_type imag);
    inline                             complex(const complex< pod_type >& c);
    
    inline       pod_type              abs() const;
    inline       void                  polar(const pod_type& rho, const pod_type& theta = 0);
    
    inline       complex< pod_type >   operator+(const complex< pod_type >& rhs);
    inline       complex< pod_type >   operator-(const complex< pod_type >& rhs);
    inline       complex< pod_type >   operator*(const complex< pod_type >& rhs);
    inline       complex< pod_type >   operator/(const complex< pod_type >& rhs);
    
    inline const complex< pod_type >&  operator+=(const complex< pod_type >& rhs);
    inline const complex< pod_type >&  operator-=(const complex< pod_type >& rhs);
    inline const complex< pod_type >&  operator*=(const complex< pod_type >& rhs);
    inline const complex< pod_type >&  operator/=(const complex< pod_type >& rhs);
    
    template< typename U > inline       complex< T >  operator+(const U& rhs);
    template< typename U > inline       complex< T >  operator-(const U& rhs);
    template< typename U > inline       complex< T >  operator*(const U& rhs);
    template< typename U > inline       complex< T >  operator/(const U& rhs);
    
    template< typename U > inline const complex< T >& operator+=(const U& rhs);
    template< typename U > inline const complex< T >& operator-=(const U& rhs);
    template< typename U > inline const complex< T >& operator*=(const U& rhs);
    template< typename U > inline const complex< T >& operator/=(const U& rhs);
};



/*!
 * @brief           Destructor for a complex number.
 * @details         Frees allocated memory
 */
template< typename T >
inline
complex< T >::~complex()
{}

/*!
 * @brief           Default constructor for creating a complex number.
 * @details         A default complex number contains only a real value
 *                  of 0 and a imaginary value of 0.
 */
template< typename T >
inline
complex< T >::complex()
: re(0)
, im(0)
{}



/*!
 * @brief           Constructs a complex number with the same real and imaginary
 *                  value.
 * @details         Constructs a complex number which has the same value for the
 *                  real and the imaginary value.
 *
 * @param[in]       real_imag The real and complex value.
 *
 */
template< typename T >
inline
complex< T >::complex(const T real_imag)
: re(real_imag)
, im(real_imag)
{}



/*!
 * @brief           Constructs a complex number with the given value for the real
 *                  value and the given value for the imaginary value.
 * @details         Constructs a complex number that has a real value of the value
 *                  that is given by the parameter _real_ and a imaginary value
 *                  given by the parameter _imag_
 *
 * @param[in]       real The real value of the complex number that is supposed to
 *                  be constructed.
 * @param[in]       imag The imaginary value of the complex number that is supposed
 *                  to be constructed.
 *
 */
template< typename T >
inline
complex< T >::complex(const T real, const T imag)
: re(real)
, im(imag)
{}


/*!
 * @brief           The copy constructor for a complex number.
 * @details         Constructs a complex number that is a copy of the given complex
 *                  number.
 *
 * @param[in]       c The complex number that is supposed to be copied.
 */
template< typename T >
inline
complex< T >::complex(const complex< T >& c)
: re(c.re)
, im(c.im)
{}

/*!
 * @brief           Returns the absolute value of the current complex number.
 * @details         The absolute value of a complex number \f$\underline{z}\f$
 *                  is defined as the euclidian norm of its real and imaginary
 *                  parts.
 *                  \f[
 *                      |\underline{z}| = \sqrt{\Re{(z)}^2 + \Im{(z)}^2}
 *                  \f]
 *
 * @return          The absolute value of the current complex number as defined
 *                  in the detailed description.
 */
template< typename T >
inline
T complex< T >::abs() const
{
    if (same_type< T, double >::value)
    {
        return sqrt(re * re + im * im);
    }
    else if (same_type< T, float >::value)
    {
        return sqrtf(re * re + im * im);
    }
    else if (same_type< T, long double >::value)
    {
        return sqrtl(re * re + im * im);
    }
    else
    {
        return sqrtf(static_cast< float >(re * re) + static_cast< float >(im * im));
    }
}



/*!
 * @brief           Constructs a complex number by polar coordinates.
 * @details         The complex number \f$\underline{z}\f$ which is defined by polar
 *                  coordinates can be expressed in real and imaginary parts by applying
 *                  \f{eqnarray*}{
 *                      \Re(\underline{z}) &=& \rho\cdot \cos(\theta)\\
 *                      \Im(\underline{z}) &=& \rho\cdot \sin(\theta)
 *                  \f}
 *                  where \f$\rho\f$ is the absolute value of the complex number in polar
 *                  coordinates and \f$\theta\f$ is the argument of the copmlex number
 *                  as a radiant angle.
 *
 * @param[in]       rho The absolute value of the complex number.
 * @param[in]       theta The argument of the complex number as radiant angle.
 */
template< typename T >
inline
void complex< T >::polar(const T& rho, const T& theta)
{
    if ( same_type< T, double >::value )
    {
        re = rho * cos(theta);
        im = rho * sin(theta);
    }
    else if ( same_type< T, float >::value )
    {
        re = rho * cosf(theta);
        im = rho * sinf(theta);
    }
    else if ( same_type< T, long double >::value )
    {
        re = rho * cosl(theta);
        im = rho * sinl(theta);
    }
    else
    {
        re = static_cast< T >(rho * cosf(static_cast< float >(theta)));
        im = static_cast< T >(rho * cosf(static_cast< float >(theta)));
    }
}

/*!
 * @brief           Addition operator for two complex numbers.
 * @details         Adds two complex nubers by adding the real values
 *                  and the imaginary values of both complex numbers.
 *
 * @param[in]       rhs Complex number on the right handside of the addition
 *                  operator
 *
 * @return          The result of the addition.
 */
template< typename T >
inline
complex< T > complex< T >::operator+(const complex< T >& rhs)
{
    return complex< T >(re + rhs.re, im + rhs.im);
}

/*!
 * @brief           Subtraction operator for two complex numbers.
 * @details         Subtracts the given complex number from the current
 *                  complex number by subtracting the real and imaginary values
 *                  of the current complex number by the real and imaginary
 *                  values of the given complex number.
 *
 * @param[in]       rhs The complex number that is supposed to be subtracted from
 *                  the current complex number.
 *
 * @return          A new complex number containing the result.
 */
template< typename T >
inline
complex< T > complex< T >::operator-(const complex< T >& rhs)
{
    return complex< T >(re - rhs.re, im - rhs.im);
}

/*!
 * @brief           Multiplication operator for two complex numbers.
 * @details         Multiplies the given complex number to the current complex
 *                  number by applying the following calculation
 *                  \f{eqnarray*}{
 *                      \underline{z} &=& \underline{x}\cdot \underline{y}\\
 *                      \underline{z} &=& (\Re(\underline{x})\cdot\Re(\underline{y})
 *                      - \Im(\underline{x})\cdot\Im(\underline{y}))+(\Im(\underline{x})\cdot\Re(\underline{y})
 *                      + \Re(\underline{x})\cdot\Im(\underline{y}))i
 *                  \f}
 *                  where \f$\underline{x},\;\underline{y},\;\underline{z}\f$
 *                  denoting complex numbers.
 *
 * @param[in]       rhs The complex number that is supposed to be multiplied to
 *                  the current complex number.
 *
 * @return          A new complex number containing the multiplication result.
 */
template< typename T >
inline
complex< T > complex< T >::operator*(const complex< T >& rhs)
{
    return complex< T >(re * rhs.re - im * rhs.im,   // real value of this number
                        im * rhs.re + re * rhs.im);  // imag value of this number
}

/*!
 * @brief           Division operator for two complex numbers.
 * @details         Divides the given complex number from the current complex
 *                  number by applying the following calculation
 *                  \f{eqnarray*}{
 *                      \underline{z} &=& \underline{x}\cdot \underline{y}\\
 *                      \underline{z} &=& \left(\frac{\Re(\underline{x})\cdot\Re(\underline{y})
 *                          - \Im(\underline{x})\cdot\Im(\underline{y})}{\Re(\underline{y})^2
 *                          + \Im(\underline{y})^2}\right) + \left(\frac{\Im(\underline{x})\cdot\Re(\underline{y})
 *                          +\Re(\underline{x})\cdot\Im(\underline{y})}{\Re(\underline{y})^2
 *                          + \Im(\underline{y})^2}\right)i
 *                  \f}
 *                  where \f$\underline{x},\;\underline{y},\;\underline{z}\f$
 *                  denoting complex numbers.
 *
 * @param[in]       rhs The complex number that is supposed to be devided by.
 *
 * @return          A new complex number containing the division result.
 */
template< typename T >
inline
complex< T > complex< T >::operator/(const complex< T >& rhs)
{
    return complex< T >((re * rhs.re + im * rhs.im) / rhs.norm(),    // real value
                        (im * rhs.re - re * rhs.im) / rhs.norm());   // imag value
}

/*!
 * @brief           Addition operator for an integer value and a complex number.
 * @details         Calculates the sum of an integer value and a complex number
 *                  by casting the int to a complex number and performing the
 *                  addition of two complex numbers.
 *
 * @param[in]       lhs The integer value on the left handside of the addition
 *                  operator.
 * @param[in]       rhs The complex value on the right handside of the addition
 *                  operator.
 *
 * @return          A new complex number containing the result of the addition.
 *
 * @sa              complex::operator+
 */
template< typename T >
complex< T > operator+(int lhs, complex< T > rhs)
{
    return complex< T >(lhs, 0) + rhs;
}

/*!
 * @brief           Addition operator for an non-complex value and a complex number.
 * @details         Calculates the sum of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  addition of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the addition
 *                  operator.
 *
 * @return          A new complex number containing the result of the addition.
 *
 * @sa              complex::operator+
 */
template< typename T >
template< typename U >
inline
complex< T > complex< T >::operator+(const U& rhs)
{
    complex< T > c(rhs, 0);
    return *this + c;
}

/*!
 * @brief           Subtraction operator for an non-complex value and a complex number.
 * @details         Calculates the difference of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  subtraction of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the subtraction
 *                  operator.
 *
 * @return          A new complex number containing the result of the subtraction.
 *
 * @sa              complex::operator-
 */
template< typename T >
template< typename U >
inline
complex< T > complex< T >::operator-(const U& rhs)
{
    complex< T > c(rhs, 0);
    return *this - c;
}

/*!
 * @brief           Multiplication operator for an non-complex value and a complex number.
 * @details         Calculates the product of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  multiplication of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the multiplication
 *                  operator.
 *
 * @return          A new complex number containing the result of the multiplication.
 *
 * @sa              complex::operator*
 */
template< typename T >
template< typename U >
inline
complex< T > complex< T >::operator*(const U& rhs)
{
    complex< T > c(rhs, 0);
    return *this * c;
}

/*!
 * @brief           Division operator for an non-complex value and a complex number.
 * @details         Calculates the quotient of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  division of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the division
 *                  operator.
 *
 * @return          A new complex number containing the result of the division.
 *
 * @sa              complex::operator/
 */
template< typename T >
template< typename U >
inline
complex< T > complex< T >::operator/(const U& rhs)
{
    complex< T > c(rhs, 0);
    return *this / c;
}



/*!
 * @brief           The addition assignment operator for two complex numbers.
 * @details         Adds a complex number to the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  addition assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to add two complex numbers
 *                  see complex::operator+
 */
template< typename T >
inline
const complex< T >& complex< T >::operator+=(const complex< T >& rhs)
{
    re += rhs.re;
    im += rhs.im;
    return *this;
}

/*!
 * @brief           The subtraction assignment operator for two complex numbers.
 * @details         Subtracts a complex number from the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  subtraction assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to subtract two complex numbers
 *                  see complex::operator-
 */
template< typename T >
inline
const complex< T >& complex< T >::operator-=(const complex< T >& rhs)
{
    re -= rhs.re;
    im -= rhs.im;
    return *this;
}

/*!
 * @brief           The multiplication assignment operator for two complex numbers.
 * @details         Multiplies a complex number to the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  multiplication assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to multiply two complex numbers
 *                  see complex::operator*
 */
template< typename T >
inline
const complex< T >& complex< T >::operator*=(const complex< T >& rhs)
{
    complex< T > c = *this * rhs;
    re = c.re;
    im = c.im;
    return *this;
}

/*!
 * @brief           The division assignment operator for two complex numbers.
 * @details         Divides a complex number by a complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  division assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to divide two complex numbers
 *                  see complex::operator/
 */
template< typename T >
inline
const complex< T >& complex< T >::operator/=(const complex< T >& rhs)
{
    complex< T > c = *this / rhs;
    re = c.re;
    im = c.im;
    return *this;
}

/*!
 * @brief           The addition assignment operator for a complex number and
 *                  a non-complex number.
 * @details         Adds a non-complex number to the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number. The
 *                  non-complex number gets casted into a complex number and
 *                  is then added to the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  addition assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to add two complex numbers
 *                  see complex::operator+
 */
template< typename T >
template< typename U >
inline
const complex< T >& complex< T >::operator+=(const U& rhs)
{
    if ( is_complex< U >::value )
    {
        *this += rhs;
    }
    else
    {
        complex< T > c(rhs, 0);
        *this += c;
    }
    return *this;
}

/*!
 * @brief           The subtraction assignment operator for a complex number and
 *                  a non-complex number.
 * @details         Subtracts a non-complex number from the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number. The
 *                  non-complex number gets casted into a complex number and
 *                  is then subtracted from the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  subtraction assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to subtract two complex numbers
 *                  see complex::operator-
 */
template< typename T >
template< typename U >
inline
const complex< T >& complex< T >::operator-=(const U& rhs)
{
    if (is_complex< U >::value == true)
    {
        *this -= rhs;
    }
    else
    {
        complex< T > c(rhs, 0);
        *this -= c;
    }
    return *this;
}

/*!
 * @brief           The mutliplication assignment operator for a complex number and
 *                  a non-complex number.
 * @details         Multiplies a non-complex number from the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number. The
 *                  non-complex number gets casted into a complex number and
 *                  is then multiplied with the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  mutliplication assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to multiply two complex numbers
 *                  see complex::operator*
 */
template< typename T >
template< typename U >
inline
const complex< T >& complex< T >::operator*=(const U& rhs)
{
    if (is_complex< U >::value == true)
    {
        *this *= rhs;
    }
    else
    {
        complex< T > c(rhs, 0);
        *this *= c;
    }
    return *this;
}

/*!
 * @brief           The division assignment operator for a complex number and
 *                  a non-complex number.
 * @details         Divides a non-complex number from the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number. The
 *                  non-complex number gets casted into a complex number and
 *                  then the division is performed.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  division assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to divide two complex numbers
 *                  see complex::operator/
 */
template< typename T >
template< typename U >
inline
const complex< T >& complex< T >::operator/=(const U& rhs)
{
    if (is_complex< U >::value == true)
    {
        *this /= rhs;
    }
    else
    {
        complex< T > c(rhs, 0);
        *this /= c;
    }
    return *this;
}


/*- Additional operators -*/


///*!
// * @brief           Multiplication operator for an non-complex value and a complex number.
// * @details         Calculates the product of an non-complex value and a complex number
// *                  by casting the non-complex to a complex number and performing the
// *                  multiplication of two complex numbers.
// *
// * @param[in]       lhs The non-complex value on the left handside of the multiplication
// *                  operator.
// * @param[in]       rhs The complex value on the right handside of the multiplication
// *                  operator.
// *
// * @return          A new complex number containing the result of the multiplication.
// *
// * @sa              complex::operator*
// */
//template<typename eT, typename T>
//complex< eT > operator*(const T& lhs, complex< eT > rhs)
//{
//    return complex< eT >(lhs, 0) * rhs;
//}
//
///*!
// * @brief           Division operator for an non-complex value and a complex number.
// * @details         Calculates the quotient of an non-complex value and a complex number
// *                  by casting the non-complex to a complex number and performing the
// *                  division of two complex numbers.
// *
// * @param[in]       lhs The non-complex value on the left handside of the division
// *                  operator.
// * @param[in]       rhs The complex value on the right handside of the division
// *                  operator.
// *
// * @return          A new complex number containing the result of the division.
// *
// * @sa              complex::operator/
// */
//template<typename eT, typename T>
//complex< eT > operator/(const T& lhs, complex< eT > rhs)
//{
//    return complex< eT >(lhs, 0) / rhs;
//}
//
///*!
// * @brief           Addition operator for an non-complex value and a complex number.
// * @details         Calculates the sum of an non-complex value and a complex number
// *                  by casting the non-complex to a complex number and performing the
// *                  addition of two complex numbers.
// *
// * @param[in]       lhs The non-complex value on the left handside of the addition
// *                  operator.
// * @param[in]       rhs The complex value on the right handside of the addition
// *                  operator.
// *
// * @return          A new complex number containing the result of the addition.
// *
// * @sa              complex::operator+
// */
//template<typename eT, typename T>
//complex< eT > operator+(const T& lhs, complex< eT > rhs)
//{
//    return complex< eT >(lhs, 0) + rhs;
//}
//
///*!
// * @brief           Subtraction operator for an non-complex value and a complex number.
// * @details         Calculates the difference of an non-complex value and a complex number
// *                  by casting the non-complex to a complex number and performing the
// *                  subtraction of two complex numbers.
// *
// * @param[in]       lhs The non-complex value on the left handside of the subtraction
// *                  operator.
// * @param[in]       rhs The complex value on the right handside of the subtraction
// *                  operator.
// *
// * @return          A new complex number containing the result of the subtraction.
// *
// * @sa              complex::operator-
// */
//template<typename eT, typename T>
//complex< eT > operator-(const T& lhs, complex< eT > rhs)
//{
//    return complex< eT >(lhs, 0) - rhs;
//}

/*!
 * @brief           Outstream operator overload to print the complex number in a nice
 *                  format to the outstream.
 * @details         The overload for the outstream operator can be used like this
 *                  @code
 *                      complex<double> z(2.5, 3.2);
 *                      std::cout << "z = " << z << std::endl;
 *
 *                      // The output looks as follows:
 *                      // z = +2.5000e+00+3.2000e+00i
 *                  @endcode
 *                  The number is printed in scientific format.
 *
 * @param[in]       o The outstream object.
 * @param[in]       c The complex number that is supposed to be printed.
 * @tparam          T The template type of the given complex number.
 *
 * @return          The reference to the outstream object.
 */
template< typename T >
std::ostream& operator<<(std::ostream& o, const complex< T >& c)
{
    std::ios::fmtflags f( std::cout.flags() );
    o << std::scientific << std::setprecision(4);
    
    o << std::noshowpos;
    o << c.re;
    o << (c.im < 0 ? " - " : " + ") << (c.im == 0 ?  0 : std::abs(c.im)) << "i";
    
    std::cout.flags( f );
    return o << std::fixed;
}

/*!
 * @}
 */

PSOFFT_END

#endif /* complex.hpp */
