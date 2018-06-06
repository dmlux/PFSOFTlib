//
//  fftw_wrapper.hpp
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

#ifndef PFSOFTlib_fftw_wrapper_hpp
#define PFSOFTlib_fftw_wrapper_hpp

PFSOFT_BEGIN

extern "C"
{
    /*- FFTW FUNCTIONS -*/
    void uzl_fftw_layer_wise_DFT2_grid3D (int cols, int rows, int lays, double* arr);
    void uzl_fftw_layer_wise_IDFT2_grid3D(int cols, int rows, int lays, double* arr);
}

PFSOFT_END

#endif /* fftw_wrapper.hpp */
