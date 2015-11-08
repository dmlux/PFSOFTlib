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

#ifndef PSOFFT_compiler_hpp
#define PSOFFT_compiler_hpp

// Library version
#undef  PSOFFT_MAJOR
#define PSOFFT_MAJOR 1

#undef  PSOFFT_MINOR
#define PSOFFT_MINOR 0

#undef  PSOFFT_PATCH
#define PSOFFT_PATCH 0

// Define the default number of used threads for multithreaded parts
#ifdef _OPENMP

    // include OpenMP headers for the OpenMP API
    #include <omp.h>

    // Number of default used threads is the max
    // number of possible threads
    #undef  PSOFFT_MAX_THREADS
    #define PSOFFT_MAX_THREADS omp_get_max_threads()
#else

    // setting the number of default used threads for
    // multithreaded parts to 1
    #undef  PSOFFT_MAX_THREADS
    #define PSOFFT_MAX_THREADS 1
#endif

/*- Processor information -*/
#define PSOFFT_PROJECT_ARCH "x86_64"

/*- SOFT config -*/
// SOFT Threshold
#undef  DSOFT_THRESHOLD
#define DSOFT_THRESHOLD 20

/*- Namespace macros -*/
// Macro shortcut for standard PSOFFT namespace
#undef  PSOFFT_BEGIN
#define PSOFFT_BEGIN           namespace psofft {

#undef  PSOFFT_END
#define PSOFFT_END             }

// Macro for nested standard namespace.
#define PSOFFT_NAMESPACE(a)    namespace psofft { namespace a {
    
#undef  PSOFFT_NAMESPACE_END
#define PSOFFT_NAMESPACE_END   }}

/*- Debugging messages -*/
// The debug mode flag is used to determine
// whether debug messages of methods should
// be printed or not. What messages get printed
// depends on the specific flags.
#undef  PSOFFT_DEBUG
#define PSOFFT_DEBUG 1

// the PSOFFT_SHOW_WARNINGS flag indicates
// that warning messages should be printed
// to console if they occure in execution.
#undef  PSOFFT_SHOW_WARNINGS
#define PSOFFT_SHOW_WARNINGS 1

// The PSOFFT_SHOW_ERRORS flag indicates
// that error messages should be printed
// to console if they occure in execution.
#undef  PSOFFT_SHOW_ERRORS
#define PSOFFT_SHOW_ERRORS 1

/*- Debugging log -*/
// Printing error message to the stderr console
#define psofft_cond_e(condition, fmt, ...)\
        do {\
            if(condition && PSOFFT_DEBUG && PSOFFT_SHOW_ERRORS) {\
                fprintf(stderr, "** [PSOFFTLib error]   %s:%d:%s(): " fmt " **\n",\
                    __FILE__, __LINE__, __func__, __VA_ARGS__\
                );\
                exit(EXIT_FAILURE);\
            }\
        } while(0)

// printing warning message to the stderr console
#define psofft_cond_w(condition, fmt, ...)\
        do {\
            if (condition && PSOFFT_DEBUG && PSOFFT_SHOW_WARNINGS)\
                fprintf(stderr, "** [PSOFFTLib warning] %s:%d:%s(): " fmt " **\n",\
                    __FILE__, __LINE__, __func__, __VA_ARGS__\
                );\
        } while(0)
        
#define psofft_cond_w_ret(condition, fmt, ...)\
        do {\
            if (condition && PSOFFT_DEBUG && PSOFFT_SHOW_WARNINGS) {\
                fprintf(stderr, "** [PSOFFTLib warning] %s:%d:%s(): " fmt " **\n",\
                    __FILE__, __LINE__, __func__, __VA_ARGS__\
                );\
                return;\
            }\
        } while(0)
      
/*- Attributes -*/
// force inlining of function
#undef  psofft_inline
#define psofft_inline          inline __attribute__((always_inline))
        
#undef  psofft_deprecated
#define psofft_deprecated      __attribute__((deprecated))

#define psofft_nonnull(...)    __attribute__((nonnull(__VA_ARGS__)))

#undef  psofft_pure
#define psofft_pure            __attribute__((pure))

#define psofft_aligned(n)      __attribute__((aligned((n))))

#endif
