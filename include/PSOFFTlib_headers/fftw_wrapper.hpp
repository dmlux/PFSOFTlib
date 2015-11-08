//
//  fftw_wrapper.hpp
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

#ifndef PSOFFTlib_fftw_wrapper_hpp
#define PSOFFTlib_fftw_wrapper_hpp

PSOFFT_BEGIN

extern "C"
{
    /*- FFTW FUNCTIONS -*/
    void uzl_fftw_layer_wise_DFT2_grid3D (int cols, int rows, int lays, double* arr);
    void uzl_fftw_layer_wise_IDFT2_grid3D(int cols, int rows, int lays, double* arr);
}

PSOFFT_END

#endif /* fftw_wrapper.hpp */
