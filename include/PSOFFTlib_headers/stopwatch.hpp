//
//  stopwatch.hpp
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

#ifndef PSOFFTlib_stop_watch_dec_hpp
#define PSOFFTlib_stop_watch_dec_hpp

PSOFFT_BEGIN

/*!
 * @brief   A measurment tool for getting execution time of code snippets.
 * @details The stop watch can be used to measure the execution time of a
 *          specified code snippet. The time can be returned in common time
 *          units like nano seconds, micro seconds, milli seconds, seconds,
 *          minutes and hours.
 *
 * @since   0.0.1
 *
 * @author  Denis-Michael Lux <denis.lux@icloud.com>
 * @date    01.05.15
 */
class stopwatch
{
    struct timeval start;   //!< The start time reference
    
    stopwatch();            
    
public:
    static stopwatch tic();
    ~stopwatch();
    
    double toc();
    double toc_micros();
    double toc_millis();
    double toc_seconds();
    double toc_minutes();
    double toc_hours();
};

PSOFFT_END

#endif /* stopwatch.hpp */
