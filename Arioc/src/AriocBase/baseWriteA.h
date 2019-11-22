/*
  baseWriteA.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#pragma once
#define __baseWriteA__

#pragma region structs
struct RowWriteCounts
{
    // unpaired
    INT64   nMapped1;
    INT64   nMapped2;
    INT64   nUnmapped;

    // paired
    INT64   nConcordantRows;
    INT64   nDiscordantRows;
    INT64   nRejectedRows;
    INT64   nUnmappedRows;
};
#pragma endregion

/// <summary>
/// Class <c>baseWriteA</c> declares common functionality for derived class tuWriteA in AriocU and AriocP
/// </summary>
class baseWriteA : public tuBaseA
{
    protected:
        QBatch*                         m_pqb;
        OutputFormatType                m_oft;
        baseARowBuilder*                m_parb;
        QAI                             m_qaiU;
        PAI                             m_paiU;
        QAI                             m_qaiU0;
        QAI                             m_qaiU1;
        INT32                           m_iBuf;
        volatile UINT32*                m_pi;
        bool                            m_wantCounts;
        AriocTaskUnitMetrics*           m_ptum;
        HiResTimer                      m_hrt;

    protected:
        baseWriteA( void );
        baseWriteA( QBatch* pqb, OutputFormatType oft, baseARowBuilder* parb, const char* tumKey, INT32 iBuf, volatile UINT32* pi, bool wantCounts );

        void main( void );

        virtual void emitRows( RowWriteCounts& rwc ) = 0;
        virtual PAI* buildTempQAI( PAI* pPAI, Qwarp* pQw ) = 0;
        virtual QAI* buildTempQAI( Qwarp* pQw, UINT32 iw, INT16 iq ) = 0;

    public:
        virtual ~baseWriteA( void );
};
