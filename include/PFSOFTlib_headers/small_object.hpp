//
//  small_ojbect.hpp
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

#ifndef PFSOFTlib_small_object_hpp
#define PFSOFTlib_small_object_hpp

PFSOFT_BEGIN

struct Chunk
{
private:
    friend class FixedAllocator;
    
    void Init(std::size_t blockSize, unsigned char blocks);
    void* Allocate(std::size_t blockSize);
    void Deallocate(void* p, std::size_t blockSize);
    void Release();
    unsigned char* pData_;
    unsigned char firstAvailableBLock_;
    unsigned char blocksAvailable_;
};

class FixedAllocator
{
    std::size_t blockSize_;
    unsigned char numBlocks_;
    typedef std::vector<Chunk> Chunks;
    Chunks chunks_;
    Chunk* allocChunk_;
    Chunk* deallocChunk_;
    
public:
    void* Allocate();
};

PFSOFT_END

#endif
