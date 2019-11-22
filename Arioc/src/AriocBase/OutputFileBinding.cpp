/*
  OutputFileBinding.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma region static member variables
WinGlobalPtr<OutputFileBinding*> OutputFileBinding::m_ofb;
AriocBase* OutputFileBinding::m_pab = NULL;
INT32 OutputFileBinding::BuffersPerBinding = 0;
INT32 OutputFileBinding::ThreadsPerGPU = 0;
INT32 OutputFileBinding::m_nBuffersPerGPU = 0;
#pragma endregion

#pragma region constructors and destructor
/// [private] default constructor
OutputFileBinding::OutputFileBinding() : m_arf(arfNone), pARW(NULL)
{
}

/// constructor
OutputFileBinding::OutputFileBinding( AriocBase* _pab, AlignmentResultFlags _arf, baseARowWriter* _parw ) : m_arf(_arf), pARW( _parw )
{
    // save a static reference to the one-and-only AriocBase instance
    m_pab = _pab;
}

/// destructor
OutputFileBinding::~OutputFileBinding()
{
}
#pragma endregion

#pragma region private methods
/// [private static] method computeOutputBufferCounts
void OutputFileBinding::computeOutputBufferCounts()
{
    /* Choose a reasonable number of threads to write alignments (see tuTailP and tuTailU):
        - Overall, there is one buffer per thread, per OutputFileBinding, per GPU.  (See OutputFileBinding::Reset().)
           The value we compute here represents only the number of buffers per OutputFileBinding instance.
        - We exclude a few threads from the computation:  one thread per GPU, plus the application's main thread.
    */

    // scan the static list for the number of OutputFileBinding instances associated with an active row writer (ARowWriter)
    INT32 nOFB = 0;
    for( size_t n=0; n<m_ofb.Count; ++n )
    {
        if( m_ofb.p[n]->pARW->IsActive )
            ++nOFB;
    }

    // sanity check
    if( nOFB == 0 )
        throw new ApplicationException( __FILE__, __LINE__, "no active output file bindings" );

    // exclude application main thread and per-GPU threads
    INT32 nExcluded = 1 + m_pab->nGPUs;

    // compute the total number of threads we are willing to allocate for writing formatted output
    INT32 nThreadsAvailable = max2( 1, m_pab->nLPs-nExcluded );

    // compute the tentative number of threads (buffers) per OutputFileBinding instance per GPU
    INT32 nThreads = nThreadsAvailable / (nOFB * m_pab->nGPUs);
    OutputFileBinding::BuffersPerBinding = min2( 10, max2( 1, nThreads ) );     // 10 threads should suffice

    // compute the tentative number of buffers per GPU
    m_nBuffersPerGPU = OutputFileBinding::BuffersPerBinding * nOFB;

    // adjust the thread count downward so as not to exceed the configured limit
    INT32 nTotalThreads = m_nBuffersPerGPU * m_pab->nGPUs;
    while( (OutputFileBinding::BuffersPerBinding > 1) && (nTotalThreads > nThreadsAvailable) )
    {
        m_nBuffersPerGPU = (--OutputFileBinding::BuffersPerBinding) * nOFB;
        nTotalThreads = m_nBuffersPerGPU * m_pab->nGPUs;
    }

    // count the number of distinct OutputFormatTypes in active OutputFileBinding instances 
    INT32 nActiveOFTs = 0;
    INT32 iOFTs = 0;
    for( size_t n=0; n<m_ofb.Count; ++n )
    {
        OutputFileBinding* pofb = m_ofb.p[n];

        if( pofb->pARW->IsActive && !(pofb->pARW->OFT & iOFTs) )
        {
            iOFTs |= pofb->pARW->OFT;
            ++nActiveOFTs;
        }
    }

    // compute the number of threads per GPU for concurrent formatting of alignment results (see tuTailP::main())
    OutputFileBinding::ThreadsPerGPU = nActiveOFTs * OutputFileBinding::BuffersPerBinding;
}
#pragma endregion

#pragma region public methods
#if TODO_CHOP_IF_UNUSED
/// [public] method CopyToOutputBuffer
void OutputFileBinding::CopyToOutputBuffer( INT32 _iBuf, WinGlobalPtr<char>& _rowBuf )
{
#pragma warning( push )
#pragma warning( disable:4996 )                 // (don't nag us about memcpy being "unsafe")

    if( this->pARW->IsActive )
    {
        char* pOutBuf = this->pARW->Reserve( m_buf.p[_iBuf], _rowBuf.n );
        memcpy( pOutBuf, _rowBuf.p, _rowBuf.n );
        this->pARW->Commit( m_buf.p[_iBuf], _rowBuf.n );
    }

#pragma warning( pop )
}
#endif

/// [public] method ReserveOutputBuffer
char* OutputFileBinding::ReserveOutputBuffer( INT32 _iBuf, UINT32 _cb )
{
    if( !this->pARW->IsActive )
        return NULL;

    return this->pARW->Reserve( m_buf.p[_iBuf], _cb );
}

/// [public] method CommitOutputBuffer
void OutputFileBinding::CommitOutputBuffer( INT32 _iBuf, UINT32 _cb )
{
    if( !this->pARW->IsActive )
        return;

    this->pARW->Commit( m_buf.p[_iBuf], _cb );
}

/// [public static] method CreateBinding
void OutputFileBinding::CreateBinding( AriocBase* _pab, AlignmentResultFlags _arf, baseARowWriter* _parw )
{
    OutputFormatType oft = _parw->OFT;
    if( OutputFileBinding::GetBinding( oft, _arf ) )
    {
        const char* sOFT = (oft == oftSAM ? "SAM" :
                            oft == oftSBF ? "SBF" :
                            oft == oftTSE ? "TSE" :
                            oft == oftKMH ? "KMH" :
                            "(unknown output format)");

        throw new ApplicationException( __FILE__, __LINE__, "redundantly defined report specification for %s file", sOFT );
    }

    // grow the static list of OutputFileBinding instances
    m_ofb.Realloc( m_ofb.Count+1, false );
    m_ofb.p[m_ofb.Count-1] = new OutputFileBinding( _pab, _arf, _parw );
}

/// [public static] method GetBinding
OutputFileBinding* OutputFileBinding::GetBinding( OutputFormatType _oft, AlignmentResultFlags _arf )
{
    /* Scan the static list of OutputFileBinding instances for the specified output format (SAM, SBF, etc.)
        and alignment result (concordant, discordant, etc.).
        
       When the specified alignment result bitmap has two or more bits set, only one matching bit suffices.
    */
    for( size_t n=0; n<m_ofb.Count; ++n )
    {
        if( (m_ofb.p[n]->m_arf & _arf) && (m_ofb.p[n]->pARW->OFT == _oft) )
            return m_ofb.p[n];
    }

    /* At this point there is no OutputFileBinding instance for the specified output format and alignment result. */

    return NULL;
}

/// [public static] method GetOutputFormatTypes
INT32 OutputFileBinding::GetOutputFormatTypes()
{
    INT32 ofts = 0;
    for( size_t n=0; n<m_ofb.Count; ++n )
        ofts |= m_ofb.p[n]->pARW->OFT;

    // return a bitmap of all OutputFormatTypes represented in the list of OutputFileBinding instances
    return ofts;
}

/// [public static] method GetOutputFormatTypeForCounts
OutputFormatType OutputFileBinding::GetOutputFormatTypeForCounts( void )
{
    /* The idea is to return an output format type that can be used to accumulate alignment-result counts.  The result
        here depends on the ordering in the OutputFormatType enumeration. */
    OutputFormatType oft = oftUnknown;

    for( size_t n=0; n<m_ofb.Count; ++n )
        oft = max2( oft, m_ofb.p[n]->pARW->OFT );

    return oft;
}

/// [public static] method InitOutputBuffers
void OutputFileBinding::InitOutputBuffers()
{
    // compute the number of threads to be used for each OutputFileBinding instance
    OutputFileBinding::computeOutputBufferCounts();

    // compute the total number of buffers across all GPUs
    INT32 nTotalBuffersPerBinding = OutputFileBinding::BuffersPerBinding * m_pab->nGPUs;

    // (re)allocate buffers for each OutputFileBinding instance
    for( size_t n=0; n<m_ofb.Count; ++n )
    {
        OutputFileBinding* pofb = m_ofb.p[n];

        // skip buffers for inactive
        if( !pofb->pARW->IsActive )
            continue;

        // if we need a new set of buffers...
        if( static_cast<INT32>(pofb->m_buf.Count) != nTotalBuffersPerBinding )
        {
            // zap existing buffer allocations for the current OutputFileBinding instance
            for( size_t iBuf=0; iBuf<pofb->m_buf.Count; ++iBuf )
                delete pofb->m_buf.p[iBuf];

            // allocate new buffers for the current OutputFileBinding instance
            pofb->m_buf.Realloc( nTotalBuffersPerBinding, false );
            for( INT32 iBuf=0; iBuf<nTotalBuffersPerBinding; ++iBuf )
                pofb->m_buf.p[iBuf] = new WinGlobalPtr<char>( OUTPUT_BUFFER_SIZE, false );
        }

        // reset the buffer offsets
        for( size_t iBuf=0; iBuf<pofb->m_buf.Count; ++iBuf )
            pofb->m_buf.p[iBuf]->n = 0;
    }
}

/// [public static] method GetARFforOFT
AlignmentResultFlags OutputFileBinding::GetARFforOFT( OutputFormatType _oft )
{
    INT32 arf = 0;

    // accumulate AlignmentResultFlags for all active (non-bitbucket) OutputFileBinding instances with the specified OutputFormatType
    for( size_t n=0; n<m_ofb.Count; ++n )
    {
        OutputFileBinding* pofb = m_ofb.p[n];

        if( pofb->pARW->IsActive && (pofb->pARW->OFT == _oft) )
            arf |= pofb->m_arf;
    }

    return static_cast<AlignmentResultFlags>(arf & arfMaskReport);
}

/// [public static] method GetBufferIndexBase
INT32 OutputFileBinding::GetBufferIndexBase( INT16 _gpuDeviceOrdinal )
{
    // return the index of the 0th buffer in an OutputFileBinding instance for the specified GPU
    return _gpuDeviceOrdinal * OutputFileBinding::BuffersPerBinding;
}

/// [public static] method Reset
void OutputFileBinding::Reset()
{
    // delete the baseARowWriter instances
    for( size_t n=0; n<m_ofb.Count; ++n )
    {
        OutputFileBinding* pofb = m_ofb.p[n];
        baseARowWriter* parw = pofb->pARW;
        if( parw )
        {
            // flush any remaining data to disk
            for( size_t ib=0; ib<pofb->m_buf.Count; ++ib )
                parw->Flush( pofb->m_buf.p[ib] );

            // null other references to the same baseARowWriter instance
            for( size_t nn=n+1; nn<m_ofb.Count; ++nn )
            {
                if( m_ofb.p[nn]->pARW == parw )
                    m_ofb.p[nn]->pARW = NULL;
            }

            // delete the baseARowWriter instance associated with the nth OutputFileBinding instance
            delete parw;
        }
    }

    // delete the OutputFileBinding instances
    for( size_t n=0; n<m_ofb.Count; ++n )
    {
        // zap allocated buffers
        OutputFileBinding* pofb = m_ofb.p[n];

        for( size_t iBuf=0; iBuf<pofb->m_buf.Count; ++iBuf )
            delete pofb->m_buf.p[iBuf];

        // zap the n'th OutputFileBinding instance
        delete pofb;
    }

    // free the list of OutputFileBinding instances
    m_ofb.Free();
}
#pragma endregion
