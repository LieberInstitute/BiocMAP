/*
  baseARowWriter.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.

  Notes:
   The idea here is that each calling thread (i.e., one per GPU) fills its own output buffer; the ARowWriter implementation
    flushes that buffer to disk in a thread-safe manner.
*/
#pragma once
#define __baseARowWriter__

// forward declaration
class AriocBase;

/// <summary>
/// Class <c>baseARowWriter</c> is the base implementation for class ARowWriter, which writes alignment data to a file.
/// This base class is instantiated explicitly as a placeholder for unused output-format types (see AriocBase::initARowWriters).
/// </summary>
class baseARowWriter
{
    private:
        static const INT32 BITBUCKET_SIZE = 10240;  // size of bit bucket referenced by LockOutputBuffer

    protected:
        AriocBase*              m_pab;
        OutputFileInfo          m_ofi;
        char                    m_filespecStub[FILENAME_MAX];
        char                    m_filespecExt[4];
        UINT32                  m_outFileSeqNo;
        RaiiFile                m_outFile;
        INT64                   m_nA;           // current number of reads in the current output file

    public:
        OutputFormatType        OFT;
        bool                    IsActive;
        INT64                   TotalA;         // total number of alignment results written to this ARowWriter object

    public:
        baseARowWriter( AriocBase* _pab, OutputFormatType _oft );
        virtual ~baseARowWriter( void );
        virtual char* Reserve( WinGlobalPtr<char>* _pbuf, UINT32 _cb );
        virtual void Commit( WinGlobalPtr<char>* _pbuf, UINT32 _cb );
        virtual void Flush( WinGlobalPtr<char>* _pbuf );
        virtual void Close( void );
};
