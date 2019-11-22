/*
  tuTailP.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma region constructor/destructor
/// [private] constructor
tuTailP::tuTailP()
{
}

/// constructor
tuTailP::tuTailP( QBatch* pqb ) : tuTail(pqb)
{
}

/// destructor
tuTailP::~tuTailP()
{
}
#pragma endregion

#pragma region virtual method implementations
/// [protected] method classifyAlignmentResults
void tuTailP::classifyAlignmentResults()
{
    CDPrint( cdpCD3, "[%d] %s...", m_pqb->pgi->deviceId, __FUNCTION__ );

#if TODO_CHOP_WHEN_DEBUGGED
    HiResTimer hrt;
#endif

    // classify the alignment results
    tuClassifyP classifyP( m_pqb );
    classifyP.Start();
    classifyP.Wait();

#if TODO_CHOP_WHEN_DEBUGGED
    if( m_pqb->pgi->deviceOrdinal == 0 )
        CDPrint( cdpCD0, "[%d] %s (batch 0x%016llx): after classifyP: %dms", m_pqb->pgi->deviceId, __FUNCTION__, m_pqb, hrt.GetElapsed( false ) );
#endif

    CDPrint( cdpCD3, "[%d] %s completed", m_pqb->pgi->deviceId, __FUNCTION__ );
}

/// [protected] method instantiateARowFormatter
baseWriteA* tuTailP::instantiateARowFormatter( OutputFormatType _oft, INT32 _iBuf, bool _wantCounts )
{
    switch( _oft )
    {
        case oftSAM:
            return new tuWriteSAM( m_pqb, _iBuf, &m_nSAM, _wantCounts );

        case oftSBF:
            return new tuWriteSBF( m_pqb, _iBuf, &m_nSBF, _wantCounts );

        case oftTSE:
            return new tuWriteTSE( m_pqb, _iBuf, &m_nTSE, _wantCounts );

        case oftKMH:
            return new tuWriteKMH( m_pqb, _iBuf, &m_nKMH, false );

        default:
            throw new ApplicationException( __FILE__, __LINE__, "unexpected OutputFormatType: %d", _oft );
    }
}
#pragma endregion
