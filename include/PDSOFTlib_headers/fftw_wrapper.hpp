//
//  fftw_wrapper.hpp
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

#ifndef PDSOFTlib_fftw_wrapper_hpp
#define PDSOFTlib_fftw_wrapper_hpp

PDSOFT_BEGIN

extern "C"
{
    /*- FFTW FUNCTIONS -*/
    void uzl_fftw_layer_wise_DFT2_grid3D(int cols, int rows, int lays, double* arr);
    void uzl_fftw_layer_wise_IDFT2_grid3D(int cols, int rows, int lays, double* arr);
}

PDSOFT_END

#endif /* fftw_wrapper.hpp */
