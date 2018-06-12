//
//  small_object.cpp
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

#include <pfsoft>
#include <cassert>
#include <vector>

PFSOFT_BEGIN

void Chunk::Init(std::size_t blockSize, unsigned char blocks)
{
    pData_ = new unsigned char[blockSize * blocks];
    firstAvailableBLock_ = 0;
    blocksAvailable_ = blocks;
    unsigned char i = 0;
    unsigned char* p = pData_;
    for (; i != blocks; p += blockSize)
    {
        *p = ++i;
    }
}

void* Chunk::Allocate(std::size_t blockSize)
{
    if (!blocksAvailable_) return 0;
    unsigned char* pResult = pData_ + (firstAvailableBLock_ * blockSize);
    // Update firstAvailableBlock_ to point to the next block
    firstAvailableBLock_ = *pResult;
    --blocksAvailable_;
    return pResult;
}

void Chunk::Deallocate(void *p, std::size_t blockSize)
{
    assert(p >= pData_);
    unsigned char* toRelease = static_cast<unsigned char*>(p);
    // Alignment check
    assert((toRelease - pData_) % blockSize == 0);
    *toRelease = firstAvailableBLock_;
    firstAvailableBLock_ = static_cast<unsigned char>((toRelease - pData_) / blockSize);
    // Truncation check
    assert(firstAvailableBLock_ == (toRelease - pData_) / blockSize);
    ++blocksAvailable_;
}

void* FixedAllocator::Allocate()
{
    if (allocChunk_ == 0 || allocChunk_->blocksAvailable_ == 0)
    {
        // No available memeory in this chuk
        // Try to find one
        Chunks::iterator i = chunks_.begin();
        for (;; ++i)
        {
            if (i == chunks_.end())
            {
                // All filled up - add a new chunk
                chunks_.reserve(chunks_.size() + 1);
                Chunk newChunk;
                newChunk.Init(blockSize_, numBlocks_);
                chunks_.push_back(newChunk);
                allocChunk_ = &chunks_.back();
                deallocChunk_ = &chunks_.back();
                break;
            }
            if (i->blocksAvailable_ > 0)
            {
                // Found a chunk
                allocChunk_ = &*i;
                break;
            }
        }
    }
    assert(allocChunk_ != 0);
    assert(allocChunk_->blocksAvailable_ > 0);
    return allocChunk_->Allocate(blockSize_);
}

PFSOFT_END
