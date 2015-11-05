//
//  fn_fourier_transforms.hpp
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

#ifndef PDSOFTlib_fn_fourier_transforms_hpp
#define PDSOFTlib_fn_fourier_transforms_hpp

PDSOFT_NAMESPACE(FourierTransforms)

// Forward fast Fourier transform on SO(3)
void DSOFT(grid3D< complex< double > > sample, DSOFTFourierCoefficients& fc, int threads = PDSOFT_MAX_THREADS);

// Inverse fast Fourier transform on SO(3)
void IDSOFT(const DSOFTFourierCoefficients& fc, grid3D< complex< double > >& synthesis, int threads = PDSOFT_MAX_THREADS);

PDSOFT_NAMESPACE_END

#endif /* fn_fourier_transforms.hpp */
