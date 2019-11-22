/*
  tuTail.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#pragma once
#define __tuTail__

/// <summary>
/// Class <c>tuTail</c> implements base functionality for classifying and writing alignment results to disk files
/// </summary>
class tuTail : public tuBaseA
{
    protected:
        QBatch*                     m_pqb;
        WinGlobalPtr<baseWriteA*>   m_writeA;   // alignment-row formatter (baseWriteA) references
        volatile UINT32             m_nSAM;     // shared counter for SAM-format threads
        volatile UINT32             m_nSBF;     // shared counter for SBF-format threads
        volatile UINT32             m_nTSE;     // shared counter for TSE-format threads
        volatile UINT32             m_nKMH;     // shared counter for KMH-format threads

    protected:
        void main( void );
        virtual void classifyAlignmentResults( void ) = 0;

    protected:
        tuTail( void );
        tuTail( QBatch* _pqb );
        virtual ~tuTail( void );

        virtual baseWriteA* instantiateARowFormatter( OutputFormatType _oft, INT32 _iBuf, bool _wantCounts ) = 0;
};
