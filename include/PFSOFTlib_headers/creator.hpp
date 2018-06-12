//
//  creator.hpp
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

#ifndef PFSOFTlib_creator_hpp
#define PFSOFTlib_creator_hpp

PFSOFT_BEGIN

template<class T>
struct OpNewCreator
{
    static T* Create()
    {
        return new T;
    }
    
protected:
    ~OpNewCreator()
    {}
};

template<class T>
struct MallocCreator
{
    static T* Create()
    {
        void* buffer = std::malloc(sizeof(T));
        if (!buffer) return 0;
        return new(buffer) T;
    }
    
protected:
    ~MallocCreator()
    {}
};

template<class T>
struct PrototypeCreator
{
    PrototypeCreator(T* pObj = 0) : pPrototype_(pObj)
    {}
    
    T* Create()
    {
        return pPrototype_ ? pPrototype_->clone() : 0;
    }
    
    T* GetPrototype()
    {
        return pPrototype_;
    }
    
    void SetPrototype(T* pObj)
    {
        pPrototype_ = pObj;
    }
    
private:
    T* pPrototype_;
    
protected:
    ~PrototypeCreator()
    {}
};

PFSOFT_END

#endif /* creator.hpp */
