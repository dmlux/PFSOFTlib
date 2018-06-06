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

ADD_EXECUTABLE(
    example
    ${PROJECT_SOURCE_DIR}/examples/example.cpp
)
TARGET_LINK_LIBRARIES( example PFSOFT )

ADD_EXECUTABLE(
    example_error
    ${PROJECT_SOURCE_DIR}/examples/example_error.cpp
)
TARGET_LINK_LIBRARIES( example_error PFSOFT )