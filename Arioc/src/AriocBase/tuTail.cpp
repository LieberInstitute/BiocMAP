/*
  tuTail.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma region constructor/destructor
/// [protected] constructor
tuTail::tuTail()
{

}

/// constructor
tuTail::tuTail( QBatch* _pqb ) : m_pqb( _pqb ), m_nSAM(0), m_nSBF(0), m_nTSE(0), m_nKMH(0)
{
}

/// destructor
tuTail::~tuTail()
{
}
#pragma endregion

#pragma region virtual method implementations
/// <summary>
/// Writes alignment results to disk files
/// </summary>
void tuTail::main()
{
    CDPrint( cdpCD3, "[%d] %s...", m_pqb->pgi->deviceId, __FUNCTION__ );

    // classify the alignment results
    classifyAlignmentResults();

    // choose one of the row-writer implementations to count read-mapping categories (mapped, unmapped)
    OutputFormatType oftCount = OutputFileBinding::GetOutputFormatTypeForCounts();

    // start one worker thread for each output buffer
    m_writeA.Realloc( OutputFileBinding::ThreadsPerGPU, false );
    m_writeA.n = 0;
    for( INT32 ioft=1; ioft<=oftLimit; ioft<<=1 )
    {
        OutputFormatType oft = static_cast<OutputFormatType>(ioft);
        if( OutputFileBinding::GetARFforOFT( oft ) )
        {
            INT32 iBuf0 = OutputFileBinding::GetBufferIndexBase( m_pqb->pgi->deviceOrdinal );
            for( INT32 i=0; i<static_cast<INT32>(OutputFileBinding::BuffersPerBinding); ++i )
            {
                // instantiate a row formatter (derived from baseWriteA) for the ith thread and buffer
                m_writeA.p[m_writeA.n] = instantiateARowFormatter( oft, iBuf0+i, (oftCount == oft) );

                // start the ith thread for the current output format type
                m_writeA.p[m_writeA.n]->Start();

                // iterate
                m_writeA.n++;
            }
        }
    }

    // wait for completion
    for( UINT32 n=0; n<m_writeA.n; ++n )
    {
        m_writeA.p[n]->Wait();
        delete m_writeA.p[n];
    }

    // release the current QBatch instance
    m_pqb->Release();

    CDPrint( cdpCD3, "[%d] %s completed", m_pqb->pgi->deviceId, __FUNCTION__ );
}
