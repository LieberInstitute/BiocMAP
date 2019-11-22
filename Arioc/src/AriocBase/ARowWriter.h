/*
  ARowWriter.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#pragma once
#define __ARowWriter__

#ifndef __baseARowWriter__
#include "baseARowWriter.h"
#endif

// forward declaration
class AriocBase;

/// <summary>
/// Class <c>ARowWriter</c> writes alignment data to a file.
/// </summary>
class ARowWriter : public baseARowWriter
{
    static const DWORD BUFFER_LOCK_TIMEOUT = 5 * 60 * 1000;     // 5 minutes


private:
    RaiiMutex m_mtx;

    public:
        ARowWriter( AriocBase* _pab, OutputFileInfo* _pofi );
        virtual ~ARowWriter( void );
        virtual char* Reserve( WinGlobalPtr<char>* _pbuf, UINT32 _cb );
        virtual void Commit( WinGlobalPtr<char>* _pbuf, UINT32 _cb );
        virtual void Flush( WinGlobalPtr<char>* _pbuf );
        virtual void Close( void );
};
