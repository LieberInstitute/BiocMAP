/*
  ARowWriter.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.

  Notes:
   The idea here is that each calling thread fills its own output buffer; the ARowWriter implementation
    flushes that buffer to disk in a thread-safe manner.
*/
#include "stdafx.h"

#pragma region constructors and destructor
/// constructor
ARowWriter::ARowWriter( AriocBase* pab, OutputFileInfo* pofi ) : baseARowWriter(pab,pofi->oft)
{
    // copy the specified output file info
    m_ofi = *pofi;

    /* Build the output filespec stub:
        - full subdirectory path from the configuration file
        - base filename (e.g. "AriocA")
        - alignment result flags
        - 3-digit file sequence number
        - filename extension

       Example:  e:/data/BulkData/SqA/yanhuang/110114/AriocA.cdru.001.sbf
    */
    strcpy_s( m_filespecStub, sizeof m_filespecStub, m_ofi.path );

    // append the base filename to the file specification stub
    RaiiDirectory::AppendTrailingPathSeparator( m_filespecStub );       // append a path separator
    strcat_s( m_filespecStub, sizeof m_filespecStub, m_ofi.baseName );  // append the base filename

    // append the alignment result flags to the file specification stub
    strcat_s( m_filespecStub, sizeof m_filespecStub, "." );

    // map bits to strings
    if( m_ofi.arf & arfReportMapped )                 // (arfReportMapped == arfReportConcordant)
    {
        if( pab->pifgQ->HasPairs )
            strcat_s( m_filespecStub, sizeof m_filespecStub, "c" );
        else
            strcat_s( m_filespecStub, sizeof m_filespecStub, "m" );
    }
    if( m_ofi.arf & arfReportUnmapped )   strcat_s( m_filespecStub, sizeof m_filespecStub, "u" );
    if( m_ofi.arf & arfReportDiscordant ) strcat_s( m_filespecStub, sizeof m_filespecStub, "d" );
    if( m_ofi.arf & arfReportRejected )   strcat_s( m_filespecStub, sizeof m_filespecStub, "r" );

    // save the filename extension for the output file specification
    switch( m_ofi.oft )
    {
        case oftSAM:
            strcpy_s( m_filespecExt, sizeof m_filespecExt, "sam" );
            break;

        case oftSBF:
        case oftTSE:
            strcpy_s( m_filespecExt, sizeof m_filespecExt, "sbf" );
            break;

        case oftKMH:
            strcpy_s( m_filespecExt, sizeof m_filespecExt, "kmh" );
            break;

        default:
            throw new ApplicationException( __FILE__, __LINE__, "unexpected output format type %d", m_ofi.oft );
    }

    // set a flag to indicate that this ARowWriter instance is active (i.e., not a bitbucket)
    this->IsActive = true;
}

/// destructor
ARowWriter::~ARowWriter()
{
    // close the file
    this->Close();
}
#pragma endregion

#pragma region public methods
/// [public] method Close
void ARowWriter::Close()
{
    m_outFile.Close();
    this->TotalA += m_nA;
    m_nA = 0;
}

/// [public] method Flush
void ARowWriter::Flush( WinGlobalPtr<char>* _pbuf )
{
    RaiiCriticalSection<ARowWriter> rcs;
    HiResTimer hrt;

    // do nothing if there is no data to write
    if( _pbuf->n == 0 )
        return;

    // if there is no currently-open output file, open a new output file
    if( m_outFile.Handle < 0 )
    {
        // increment the file sequence number
        ++m_outFileSeqNo;

        // build a complete file specification
        char filespec[FILENAME_MAX];
        sprintf_s( filespec, sizeof filespec, "%s.%03d.%s", m_filespecStub, m_outFileSeqNo, m_filespecExt );

        // open and initialize the file
        m_outFile.Open( filespec );
        if( m_ofi.oft == oftSAM )
            m_pab->SAMhdb.WriteHeader( m_outFile );
    }

    // write the data
    INT64 cbWritten = m_outFile.Write( _pbuf->p, _pbuf->n );
    if( cbWritten != _pbuf->n )
        throw new ApplicationException( __FILE__, __LINE__, "write failed for %s: %lld/%u bytes written", m_outFile.FileSpec.p, cbWritten, _pbuf->n );

    // reset the buffer byte counter
    _pbuf->n = 0;

    // close the current output file if we have reached the user-specified per-file limit for the number of alignments to write
    if( m_nA >= m_ofi.maxA )
        this->Close();

    // performance metrics
    AriocBase::aam.ms.WriteA += hrt.GetElapsed( false );        // (we're in a critical section so interlocked add is unnecessary)
}

/// [public] method Reserve
char* ARowWriter::Reserve( WinGlobalPtr<char>* _pbuf, UINT32 _cb )
{
    /* Each calling thread has its own buffer, so no thread synchronization here. */

    // if the buffer capacity is exceeded...
    if( (_pbuf->n + _cb) > _pbuf->cb )
    {
        // flush the buffer to disk and reset m_outBuf[iBuf].n
        this->Flush( _pbuf );
    }

    /* At this point we can use cb bytes starting at offset m_outBuf[iBuf].n in the buffer. */
    return _pbuf->p + _pbuf->n;
}

/// [public] method Commit
void ARowWriter::Commit( WinGlobalPtr<char>* _pbuf, UINT32 _cb )
{
    /* Each calling thread has its own buffer, so no thread synchronization here. */

    // count the number of bytes committed
    _pbuf->n += _cb;

    // count the number of commits (i.e., the number of alignments written to the buffer)
    ++m_nA;
}
#pragma endregion
