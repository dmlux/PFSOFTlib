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

template <typename T, typename U>
inline
matrix<T> convert(const matrix<U>& m)
{
    matrix<T> result(m.rows, m.cols);
    
    assert(m.rows == result.rows);
    assert(m.cols == result.cols);
    
    for (size_t i = 0; i < result.rows; i++)
    {
        for (size_t j = 0; j < result.cols; j++)
        {
            result(i, j) = m(i, j);
        }
    }
    
    return result;
}

template <typename T, typename U>
inline
vector<complex<T> > convert(const vector<complex<U> >& v)
{
    vector<complex<T> > result(v.size);
    
    assert(result.size == v.size);
    
    if (v.type == vector<complex<U> >::ROW)
    {
        access::rw(result.type) = vector<complex<T> >::ROW;
    }
    else
    {
        access::rw(result.type) = vector<complex<T> >::COLUMN;
    }
    
    for (size_t i = 0; i < v.size; i++)
    {
        access::rw(result[i].re) = static_cast<T>(v[i].re);
        access::rw(result[i].im) = static_cast<T>(v[i].im);
    }
    
    return result;
}

// Forward fast Fourier transform on SO(3)
void DSOFT(grid3D< complex< double > > sample, DSOFTFourierCoefficients& fc, int threads = PFSOFT_MAX_THREADS);

// Inverse fast Fourier transform on SO(3)
void IDSOFT(const DSOFTFourierCoefficients& fc, grid3D< complex< double > >& synthesis, int threads = PFSOFT_MAX_THREADS);

PFSOFT_NAMESPACE_END

#endif /* fn_fourier_transforms.hpp */
