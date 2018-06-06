//
//  fn_fourier_transforms.hpp
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

#ifndef PFSOFTlib_fn_fourier_transforms_hpp
#define PFSOFTlib_fn_fourier_transforms_hpp

PFSOFT_NAMESPACE(FourierTransforms)

// Forward fast Fourier transform on SO(3)
void DSOFT(grid3D< complex< double > > sample, DSOFTFourierCoefficients& fc, int threads = PFSOFT_MAX_THREADS);

// Inverse fast Fourier transform on SO(3)
void IDSOFT(const DSOFTFourierCoefficients& fc, grid3D< complex< double > >& synthesis, int threads = PFSOFT_MAX_THREADS);

PFSOFT_NAMESPACE_END

#endif /* fn_fourier_transforms.hpp */
