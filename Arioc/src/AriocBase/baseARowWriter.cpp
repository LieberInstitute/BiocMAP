/*
  baseARowWriter.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma region constructors and destructor
/// constructor
baseARowWriter::baseARowWriter( AriocBase* pab, OutputFormatType _oft ) : m_pab(pab),
                                                                          m_ofi("","","",oftUnknown,arfNone,0),
                                                                          m_filespecStub(), m_filespecExt(),   // (parentheses cause zero-initialization)
                                                                          m_outFileSeqNo(0), m_nA(0),
                                                                          OFT(_oft),
                                                                          IsActive(false), TotalA(0)
{
}

/// destructor
baseARowWriter::~baseARowWriter()
{
}
#pragma endregion

#pragma region virtual base method implementations
/// [public] method Reserve
char* baseARowWriter::Reserve( WinGlobalPtr<char>* _pbuf, UINT32 cb )
{
    static char bitBucket[BITBUCKET_SIZE];

    return bitBucket;
}

/// [public] method Commit
void baseARowWriter::Commit( WinGlobalPtr<char>* _pbuf, UINT32 cb )
{
    // count the total number of reported alignments
    InterlockedIncrement( reinterpret_cast<volatile UINT64*>(&this->TotalA) );
}

/// [public] method Flush
void baseARowWriter::Flush( WinGlobalPtr<char>* _pbuf )
{
}

/// [public] method Close
void baseARowWriter::Close()
{
}
#pragma endregion
