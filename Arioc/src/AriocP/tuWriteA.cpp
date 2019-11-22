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
/// Writes paired-end alignment results.
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
/// [private] method buildTempQAI (unpaired)
QAI* tuWriteA::buildTempQAI( Qwarp* pQw, UINT32 iw, INT16 iq )
{
    throw new ApplicationException( __FILE__, __LINE__, "not implemented" );
}

/// [private] method buildTempQAI
PAI* tuWriteA::buildTempQAI( PAI* pPAI, Qwarp* pQw )
{
    // make a copy of the pair alignment info
    m_paiU = *pPAI;

    // create temporary QAI instances that can be used to reference the Q sequence data and metadata
    if( pPAI->pQAI1 == NULL )
    {
        m_paiU.pQAI1 = &m_qaiU0;
        m_qaiU0.N = pQw->N[m_paiU.iq];
        m_qaiU0.qid = PACK_QID(m_paiU.iw, m_paiU.iq);
        m_qaiU0.flags = static_cast<QAIflags>(((pPAI->arf & arfPrimary) ? qaiBest : qaiNone) | ((pQw->sqId[pPAI->iq] & AriocDS::SqId::MaskMateId) ? qaiParity : qaiNone));
    }

    if( pPAI->pQAI2 == NULL )
    {
        m_paiU.pQAI2 = &m_qaiU1;
        m_qaiU1.N = pQw->N[m_paiU.iq+1];
        m_qaiU1.qid = PACK_QID( m_paiU.iw, m_paiU.iq+1 );
        m_qaiU1.flags = static_cast<QAIflags>(((pPAI->arf & arfPrimary) ? qaiBest : qaiNone) | ((pQw->sqId[pPAI->iq+1] & AriocDS::SqId::MaskMateId) ? qaiParity : qaiNone));
    }

    return &m_paiU;
}

/// [private] method emitRows
void tuWriteA::emitRows( RowWriteCounts& rwc )
{
    OutputFileBinding* pofbConcordant = OutputFileBinding::GetBinding( m_oft, arfConcordant );
    OutputFileBinding* pofbDiscordant = OutputFileBinding::GetBinding( m_oft, arfDiscordant );
    OutputFileBinding* pofbRejected = OutputFileBinding::GetBinding( m_oft, arfRejected );
    OutputFileBinding* pofbUnmapped = OutputFileBinding::GetBinding( m_oft, arfUnmapped );

    // traverse the list of paired-end reads
    UINT32 n = InterlockedExchangeAdd( m_pi, 1 );
    while( n < m_pqb->ApBuffer.n )
    {
        PAI* pPAI = m_pqb->ApBuffer.p + n;
        if( pPAI->arf & (arfWriteMate1|arfWriteMate2) )
        {
            // point to the Qwarp that contains the pair of reads
            Qwarp* pQw = m_pqb->QwBuffer.p + pPAI->iw;
            INT64 sqId = pQw->sqId[pPAI->iq];

            // write the pair
            switch( pPAI->arf & arfMaskReport )
            {
                case arfConcordant:     // two mapped mates representing a concordant paired-end mapping
                    rwc.nConcordantRows += m_parb->WriteRowPc( pofbConcordant, m_iBuf, sqId, pPAI );
                    break;

                case arfDiscordant:     // two mapped mates representing a discordant paired-end mapping
                    rwc.nDiscordantRows += m_parb->WriteRowPd( pofbDiscordant, m_iBuf, sqId, pPAI );
                    break;

                case arfRejected:       // two mapped mates, but not a concordant or discordant mapping
                    rwc.nRejectedRows += m_parb->WriteRowPr( pofbRejected, m_iBuf, sqId, buildTempQAI( pPAI, pQw ) );
                    break;

                case arfUnmapped:       // one or two unmapped mates
                    rwc.nUnmappedRows += m_parb->WriteRowPu( pofbUnmapped, m_iBuf, sqId, buildTempQAI( pPAI, pQw ) );
                    break;

                default:
                    throw new ApplicationException( __FILE__, __LINE__, "unexpected alignment reporting flags = 0x%02x", pPAI->arf & arfMaskReport );
            }
        }

        // iterate
        n = InterlockedExchangeAdd( m_pi, 1 );
    }

    // performance metrics
    if( m_wantCounts )
    {
        InterlockedExchangeAdd( &AriocBase::aam.sam.p.nRowsConcordant, rwc.nConcordantRows );
        InterlockedExchangeAdd( &AriocBase::aam.sam.p.nRowsDiscordant, rwc.nDiscordantRows );
        InterlockedExchangeAdd( &AriocBase::aam.sam.p.nRowsRejected, rwc.nRejectedRows );
        InterlockedExchangeAdd( &AriocBase::aam.sam.p.nRowsUnmapped, rwc.nUnmappedRows );
    }
}
#pragma endregion
