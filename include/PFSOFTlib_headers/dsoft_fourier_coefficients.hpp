//
//  dsoft_fourier_coefficients.hpp
//  PFSOFTlib
//
//   Created by Denis-Michael Lux on 05. November 2015.
//
//   This file is part of PFSOFTlib.
//
//   PFSOFTlib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   PFSOFTlib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with PFSOFTlib.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef PFSOFTlib_dsoft_fourier_coefficients_hpp
#define PFSOFTlib_dsoft_fourier_coefficients_hpp

PFSOFT_BEGIN

/*!
 * @brief       Collection of functions and classes for DSOFT Fourier coefficients
 *              container
 * @defgroup    DSOFTFourierCoefficients DSOFT Fourier coefficients container
 * @{
 */

/*!
 * @brief       A datastructure to manage fourier coefficients that
 *              where produced by the DSOFT algorithm described by
 *              P. J. Kostelec and D. N. Rockmore in the paper
 *              'FFTs on the Rotation Group'
 * @details     This class provides a manager class that organizes
 *              and stores the fourier coefficients that where produced
 *              by the DSOFT algorithm. The fourier coefficients are
 *              indexed over three parameter.
 *
 * @since       0.0.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        05.05.15
 */
struct DSOFTFourierCoefficients
{
private:
    matrix< complex< double > >* mem; //!< Coefficients storage
    
public:
    // public ivars
    const int bandwidth;              //!< Bandwidth of function
    
    // constructors
    DSOFTFourierCoefficients();
    DSOFTFourierCoefficients(int bandlimit);
    
    // destructor
    ~DSOFTFourierCoefficients();
    
    // methods
          complex< double >& operator()(const int& l, const int& M, const int& Mp);
    const complex< double >& operator()(const int& l, const int& M, const int& Mp) const;
    
    // prototype for the overloaded stream operator
    friend std::ostream& operator<<(std::ostream& o, const DSOFTFourierCoefficients& fc);
};

/*!
 * @}
 */

PFSOFT_END

#endif /* dsoft_fourier_coefficients.hpp */
