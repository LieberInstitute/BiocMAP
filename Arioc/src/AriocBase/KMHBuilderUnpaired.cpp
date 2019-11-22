/*
  KMHBuilderUnpaired.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma warning ( push )
#pragma warning( disable:4996 )     // (don't nag us about strcpy being "unsafe")

#pragma region constructors and destructor
/// <summary>
/// Reports kmer-hashed unpaired read sequences in SQL Server binary format.
/// </summary>
KMHBuilderUnpaired::KMHBuilderUnpaired( QBatch* pqb ) : KMHBuilderBase(pqb)
{
}

/// [public] destructor
KMHBuilderUnpaired::~KMHBuilderUnpaired()
{
}
#pragma endregion

#pragma region virtual method implementations
/// [protected] method emitBRLEQ
UINT32 KMHBuilderUnpaired::emitBRLEQ( char* pbuf, INT64 sqId, QAI* pQAI )
{
    return emitBRLEQfr( pbuf, pQAI, 0, pQAI->qid );
}

/// [protected] method emitRGID
UINT32 KMHBuilderUnpaired::emitRGID( char* pbuf, INT64 sqId, QAI* pQAI )
{
    UINT8 rgId = getRGIDfromReadMetadata( pQAI, 0, pQAI->qid );
    return intToSBF<UINT8>( pbuf, rgId );
}

/// [public] method WriteRowUm
INT64 KMHBuilderUnpaired::WriteRowUm( OutputFileBinding* pofb, INT32 iBuf, INT64 sqId, QAI* pQAI )
{
    char* p0 = pofb->ReserveOutputBuffer( iBuf, m_cbRow );
    if( p0 )
    {
        sqId = setSqIdSecBit( sqId, pQAI );

        // point to the output buffer
        char* p = p0;

        for( INT16 f=0; f<m_nEmitFields; ++f )
            p += (this->*(m_emitField[f]))(p, sqId, pQAI);     // emit the f'th field

        // compute the total number of bytes required for the row
        UINT32 cb = static_cast<UINT32>(p - p0);

        // copy the rows to the output buffer
        pofb->CommitOutputBuffer( iBuf, cb );
    }

    // return the number of rows
    return 1;
}

/// [public] method WriteRowUu
INT64 KMHBuilderUnpaired::WriteRowUu( OutputFileBinding* pofb, INT32 iBuf, INT64 sqId, QAI* pQAI )
{
    char* p0 = pofb->ReserveOutputBuffer( iBuf, m_cbRow );
    if( p0 )
    {
        // point to the output buffer
        char* p = p0;

        for( INT16 f = 0; f<m_nEmitFields; ++f )
            p += (this->*(m_emitField[f]))(p, sqId, pQAI);     // emit the f'th field

        // compute the total number of bytes required for the row
        UINT32 cb = static_cast<UINT32>(p - p0);

        // copy the rows to the output buffer
        pofb->CommitOutputBuffer( iBuf, cb );
    }

    // return the number of rows
    return 1;
}
#pragma warning ( pop )
#pragma endregion
