/*
  tuWriteA.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#pragma once
#define __tuWriteA__

/// <summary>
/// Class <c>tuWriteA</c> writes paired alignment results
/// </summary>
class tuWriteA : public baseWriteA
{
    private:
        tuWriteA( void );
        virtual void emitRows( RowWriteCounts& rwc );
        virtual QAI* buildTempQAI( Qwarp* pQw, UINT32 iw, INT16 iq );
        virtual PAI* buildTempQAI( PAI* pPAI, Qwarp* pQw );

    public:
        tuWriteA( QBatch* pqb, OutputFormatType oft, baseARowBuilder* parb, const char* tumKey, INT32 iBuf, volatile UINT32* pi, bool wantCounts );
        virtual ~tuWriteA( void );
};

/// <summary>
/// Class <c>tuWriteSAM</c> writes paired alignment results to SAM-formatted files
/// </summary>
class tuWriteSAM : public tuWriteA
{
private:
    SAMBuilderPaired    m_sbp;

public:
    tuWriteSAM( QBatch* pqb, INT32 iBuf, volatile UINT32* pi, bool wantCounts ) : 
                    tuWriteA( pqb, oftSAM, &m_sbp, "tuWriteSAM", iBuf, pi, wantCounts ), m_sbp( pqb )
    {
    }

    virtual ~tuWriteSAM( void )
    {
    }
};

/// <summary>
/// Class <c>tuWriteSBF</c> writes paired alignment results to SBF-formatted files (SAM-like fields)
/// </summary>
class tuWriteSBF : public tuWriteA
{
private:
    SBFBuilderPaired    m_sbp;

public:
    tuWriteSBF( QBatch* pqb, INT32 iBuf, volatile UINT32* pi, bool wantCounts ) :
                    tuWriteA( pqb, oftSBF, &m_sbp, "tuWriteSBF", iBuf, pi, wantCounts ), m_sbp( pqb )
    {
    }

    virtual ~tuWriteSBF( void )
    {
    }
};

/// <summary>
/// Class <c>tuWriteTSE</c> writes paired alignment results to SBF-formatted files (Terabase Search Engine fields)
/// </summary>
class tuWriteTSE : public tuWriteA
{
private:
    TSEBuilderPaired    m_tbp;

public:
    tuWriteTSE( QBatch* pqb, INT32 iBuf, volatile UINT32* pi, bool wantCounts ) :
                    tuWriteA( pqb, oftTSE, &m_tbp, "tuWriteTSE", iBuf, pi, wantCounts ), m_tbp( pqb )
    {
    }

    virtual ~tuWriteTSE( void )
    {
    }
};

/// <summary>
/// Class <c>tuWriteKMH</c> writes kmer-hashed paired read sequences
/// </summary>
class tuWriteKMH : public tuWriteA
{
private:
    KMHBuilderPaired    m_kbp;

public:
    tuWriteKMH( QBatch* pqb, INT32 iBuf, volatile UINT32* pi, bool wantCounts ) :
                    tuWriteA( pqb, oftKMH, &m_kbp, "tuWriteKMH", iBuf, pi, wantCounts ), m_kbp( pqb )
    {
    }

    virtual ~tuWriteKMH( void )
    {
    }
};

