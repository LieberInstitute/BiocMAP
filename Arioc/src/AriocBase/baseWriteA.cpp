/*
  baseWriteA.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma region constructor/destructor
/// [private] constructor
baseWriteA::baseWriteA()
{
}

/// <summary>
/// Writes alignment results.
/// </summary>
baseWriteA::baseWriteA( QBatch* pqb, OutputFormatType oft, baseARowBuilder* parb, const char* tumKey, INT32 iBuf, volatile UINT32* pi, bool wantCounts ) :
                            m_pqb( pqb ), m_oft(oft),
                            m_parb( parb ),
                            m_qaiU(), m_paiU(), m_qaiU0(), m_qaiU1(),   // (parentheses cause zero initialization)
                            m_iBuf(iBuf), m_pi(pi), m_wantCounts( wantCounts ),
                            m_ptum( AriocBase::GetTaskUnitMetrics( tumKey ) )
{
}

/// destructor
baseWriteA::~baseWriteA()
{
}
#pragma endregion

#pragma region protected methods
/// <summary>
/// Write alignment results
/// </summary>
void baseWriteA::main()
{

    CDPrint( cdpCD3, "[%d] %s (%s)...", m_pqb->pgi->deviceId, __FUNCTION__, m_ptum->Key );

    if( m_oft & OutputFileBinding::GetOutputFormatTypes() )
    {
        RowWriteCounts rwc = { 0 };
        emitRows( rwc );

        // performance metrics
        InterlockedExchangeAdd( &m_ptum->ms.Elapsed, m_hrt.GetElapsed( false ) );
    }

    CDPrint( cdpCD3, "[%d] %s (%s) completed", m_pqb->pgi->deviceId, __FUNCTION__, m_ptum->Key );
}
#pragma endregion

