//
//  compile_time_assertion.hpp
//  PFSOFTlib
//
//   Created by Denis-Michael Lux on 09. June 2018.
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

#ifndef PFSOFTlib_compile_time_assertion_hpp
#define PFSOFTlib_compile_time_assertion_hpp

PFSOFT_BEGIN

template<int> struct CompileTimeError;
template<> struct CompileTimeError<true> {};

PFSOFT_END

#define PFSOFT_STATIC_CHECK(expr, msg) \
    { pfsoft::CompileTimeError<((expr) != 0)> Error_##msg; (void)Error_##msg; }

#endif
