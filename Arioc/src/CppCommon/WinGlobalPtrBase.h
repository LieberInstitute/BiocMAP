/*
   WinGlobalPtrBase.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#pragma once
#define __WinGlobalPtrBase__

/// <summary>
/// Class <c>WinGlobalPtrBase</c> defines functionality common to all specializations of WinGlobalPtr<T>
/// </summary>
class WinGlobalPtrBase
{
protected:
    bool            m_initZero; // flag set when the allocated memory is to be zeroed
    bool            m_dtor;     // flag set when a WinGlobalPtr instance is executing its destructor

public:
    size_t          cb;         // number of allocated bytes
    size_t          Count;      // number of elements of sizeof(T)
    volatile UINT32 n;          // (unused by the WinGlobalPtr implementation)

protected:
    WinGlobalPtrBase() : m_initZero(false), m_dtor(false), cb(0), Count(0), n(0)
    {
    }

    WinGlobalPtrBase( const size_t nElements, const bool zeroAllocatedMemory ) : m_initZero(zeroAllocatedMemory), m_dtor(false), cb(0), Count(nElements), n(0)
    {
    }

public:
    virtual ~WinGlobalPtrBase( void )
    {
    }

    virtual void Realloc( const size_t, const bool ) = 0;
    virtual void Free( void ) = 0;
    virtual void New( const size_t, const bool ) = 0;
    virtual const char* GetCudaAllocationTypeName( void ) = 0;
    virtual void* GetBufferPointer( void ) = 0;
};
