/*
  tuWriteA.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma region constructor/destructor
/// [private] constructor
tuWriteA::tuWriteA()
{
}

/// <summary>
/// Writes unpaired alignment results.
/// </summary>
tuWriteA::tuWriteA( QBatch* pqb, OutputFormatType oft, baseARowBuilder* parb, const char* tumKey, INT32 iBuf, volatile UINT32* pi, bool wantCounts ) :
                        baseWriteA( pqb, oft, parb, tumKey, iBuf, pi, wantCounts )
{
}

/// destructor
tuWriteA::~tuWriteA()
{
}
#pragma endregion

#pragma region virtual method implementations
/// [private] method buildTempQAI
QAI* tuWriteA::buildTempQAI( Qwarp* pQw, UINT32 iw, INT16 iq )
{
    // create a temporary QAI instance that can be used to reference the Q sequence data and metadata
    m_qaiU.N = pQw->N[iq];
    m_qaiU.qid = PACK_QID( iw, iq );

    return &m_qaiU;
}

/// [private] method buildTempQAI (paired)
PAI* tuWriteA::buildTempQAI( PAI* pPAI, Qwarp* pQw )
{
    throw new ApplicationException( __FILE__, __LINE__, "not implemented" );
}

/// [private] method emitRows
void tuWriteA::emitRows( RowWriteCounts& rwc )
{
    OutputFileBinding* pofbMapped = OutputFileBinding::GetBinding( m_oft, arfMapped );
    OutputFileBinding* pofbUnmapped = OutputFileBinding::GetBinding( m_oft, arfUnmapped );

    // traverse the list of Qwarps in the current batch
    Qwarp* pQw = m_pqb->QwBuffer.p;
    for( UINT32 iw=0; iw<m_pqb->QwBuffer.n; ++iw )
    {
        // traverse the list of Q sequences in the Qwarp
        for( INT16 iq=0; iq<pQw->nQ; ++iq )
        {
            // write the mappings for the iq'th Q sequence in the Qwarp
            UINT32 qid = AriocDS::QID::Pack(iw,iq);
            QAIReference* pqair = m_pqb->ArBuffer.p + qid;

            INT64 sqId = pQw->sqId[iq];
            QAI* pQAI = m_pqb->AaBuffer.p + pqair->ofsQAI;
            if( pqair->nAq )
            {
                // write the mappings for the iq'th Q sequence
                for( UINT32 n=0; n<pqair->nAq; ++n )
                {
                    m_parb->WriteRowUm( pofbMapped, m_iBuf, sqId, pQAI );
                    pQAI++ ;
                }

                if( pqair->nAq > 1 )
                    rwc.nMapped2 += pqair->nAq;
                else
                    rwc.nMapped1++ ;
            }
            else
            {
                buildTempQAI( pQw, iw, iq );
                rwc.nUnmapped += m_parb->WriteRowUu( pofbUnmapped, m_iBuf, sqId, &m_qaiU );
            }
        }

        // point to the next Qwarp
        pQw++ ;
    }

    // performance metrics
    if( m_wantCounts )
    {
        InterlockedExchangeAdd( &AriocBase::aam.sam.u.nReads_Mapped1, rwc.nMapped1 );
        InterlockedExchangeAdd( &AriocBase::aam.sam.u.nReads_Mapped2, rwc.nMapped2 );
        InterlockedExchangeAdd( &AriocBase::aam.sam.u.nReads_Unmapped, rwc.nUnmapped );
    }
}
#pragma endregion
