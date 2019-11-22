/*
  tuTailU.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma region constructor/destructor
/// [private] constructor
tuTailU::tuTailU()
{
}

/// constructor
tuTailU::tuTailU( QBatch* pqb ) : tuTail(pqb)
{
}

/// destructor
tuTailU::~tuTailU()
{
}
#pragma endregion

#pragma region virtual method implementations
/// [protected] method classifyAlignmentResults
void tuTailU::classifyAlignmentResults()
{
    CDPrint( cdpCD3, "[%d] %s...", m_pqb->pgi->deviceId, __FUNCTION__ );

    // classify the alignment results
    tuClassifyU classifyU( m_pqb );
    classifyU.Start();
    classifyU.Wait();

    CDPrint( cdpCD3, "[%d] %s completed", m_pqb->pgi->deviceId, __FUNCTION__ );
}

/// [protected] method instantiateARowFormatter
baseWriteA* tuTailU::instantiateARowFormatter( OutputFormatType _oft, INT32 _iBuf, bool _wantCounts )
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
