#   Created by Denis-Michael Lux on 05. November 2015.
#
#   This file is part of PFSOFTlib.
#
#   PFSOFTlib is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   PFSOFTlib is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with PFSOFT.  If not, see <http://www.gnu.org/licenses/>.

# Creating binary and linking library
ADD_EXECUTABLE(         benchmark_dsoft_for ${PROJECT_SOURCE_DIR}/benchmark/benchmark_dsoft_for.cpp                           )
TARGET_LINK_LIBRARIES(  benchmark_dsoft_for PFSOFT                                                                            )

ADD_EXECUTABLE(         benchmark_dsoft_for_speedup ${PROJECT_SOURCE_DIR}/benchmark/benchmark_dsoft_for_speedup.cpp           )
TARGET_LINK_LIBRARIES(  benchmark_dsoft_for_speedup PFSOFT                                                                    )

ADD_EXECUTABLE(         benchmark_dsoft_for_efficiency ${PROJECT_SOURCE_DIR}/benchmark/benchmark_dsoft_for_efficiency.cpp     )
TARGET_LINK_LIBRARIES(  benchmark_dsoft_for_efficiency PFSOFT                                                                 )

ADD_EXECUTABLE(         benchmark_dsoft_inv ${PROJECT_SOURCE_DIR}/benchmark/benchmark_dsoft_inv.cpp                           ) 
TARGET_LINK_LIBRARIES(  benchmark_dsoft_inv PFSOFT                                                                            )

ADD_EXECUTABLE(         benchmark_dsoft_inv_speedup ${PROJECT_SOURCE_DIR}/benchmark/benchmark_dsoft_inv_speedup.cpp           )
TARGET_LINK_LIBRARIES(  benchmark_dsoft_inv_speedup PFSOFT                                                                    )

ADD_EXECUTABLE(         benchmark_dsoft_inv_efficiency ${PROJECT_SOURCE_DIR}/benchmark/benchmark_dsoft_inv_efficiency.cpp     )
TARGET_LINK_LIBRARIES(  benchmark_dsoft_inv_efficiency PFSOFT                                                                 )

ADD_EXECUTABLE(         benchmark_dwt_accuracy ${PROJECT_SOURCE_DIR}/benchmark/benchmark_dwt_accuracy.cpp                     )
TARGET_LINK_LIBRARIES(  benchmark_dwt_accuracy PFSOFT                                                                         )

ADD_EXECUTABLE(         benchmark_dsoft_runtime_accuracy ${PROJECT_SOURCE_DIR}/benchmark/benchmark_dsoft_runtime_accuracy.cpp )
TARGET_LINK_LIBRARIES(  benchmark_dsoft_runtime_accuracy PFSOFT                                                               )
