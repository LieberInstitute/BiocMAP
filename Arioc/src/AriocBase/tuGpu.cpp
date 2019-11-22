/*
  tuGpu.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma region constructor/destructor
/// [private] default constructor
tuGpu::tuGpu()
{
}

/// <param name="iPart">a 0-based ordinal (assigned by the caller) (one per GPU)</param>
/// <param name="pab">a reference to the application's <c>AriocP</c> or <c>AriocU</c> instance</param>
tuGpu::tuGpu( INT16 gpuDeviceOrdinal, AriocBase* pab ) : m_gpuDeviceOrdinal(gpuDeviceOrdinal),
                                                         m_pab(pab),
                                                         m_cbCgaReserved(CUDAMINRESERVEDGLOBALMEMORY)
{
    // look for the optional "cgaReserved" Xparam
    INT32 i = pab->paamb->Xparam.IndexOf( "cgaReserved" );
    if( i >= 0 )
        m_cbCgaReserved = max2( CUDAMINRESERVEDGLOBALMEMORY, pab->paamb->Xparam.Value(i) );
}

/// destructor
tuGpu::~tuGpu()
{
}
#pragma endregion

#pragma region private methods
/// [private] method sniffJeol
size_t tuGpu::sniffJeol( CudaGmemLayout* _pcgl, UINT32* _pHostBuf )
{
    /* The idea here is to ensure that the last J list in the partition for the current GPU device (identified
        by m_gpuDeviceOrdinal) is complete.  At the point where this method is called, the partition size is only
        an estimate computed by dividing the J-table size by the number of GPUs, i.e., the number of partitions.
        (See AriocBase::initHJBufferLayout.)

       Here we ensure that the last J list in the partition is not truncated by growing the partition as needed
        to fully contain the last J list.  This introduces a bit of overlap with the subsequent partition.

       The J table layout looks like this, where each 5-byte value in the table is either a J value or (if
        a list contains more than HVALUE_MAX_NJ J values) a list count:


               H table                 J table
            +-----------+          +------------+
            |    nJ=2   | -------> | Jvalue[0]  |
            +-----------+          +------------+
            |    nJ=0   | ----+    | Jvalue[1]  |
            +-----------+     |    +------------+
            |    nJ=3   | -+  +--> |     50     |
            +-----------+  |       +------------+
                  .        |       | Jvalue[0]  |
                  .        |       +------------+
                  .        .             .
                           .             .
                           .             .
                           |       +------------+
                           |       | Jvalue[49] |
                           |       +------------+
                           +-----> | Jvalue[0]  |
                                   +------------+
                                   | Jvalue[1]  |
                                   +------------+
                                   | Jvalue[2]  |
                                   +------------+
                                         .
                                         .
                                         .

       To ensure that the final J list in the buffer is not truncated, we rely on the ascending sort order of the J values
        in each list.
    */
    CudaGmemLayout::CudaGmemPartitionInfo* pcgpi = _pcgl->part + m_gpuDeviceOrdinal;

    // maximum J-list count
    const UINT64 maxnJ = (reinterpret_cast<JtableHeader*>(_pHostBuf))->maxnJ;

    // a handy mask for extracting a Jvalue5 from an 8-byte UINT64
    const UINT64 j5mask = (static_cast<UINT64>(Jvalue5::bfMaxVal_subId) << (Jvalue5::bfSize_J+Jvalue5::bfSize_s)) |
                          (static_cast<UINT64>(Jvalue5::bfMaxVal_s)     << Jvalue5::bfSize_J) |
                          static_cast<UINT64>(Jvalue5::bfMaxVal_J);

    // compute the range of 5-byte offsets into the J table buffer
    size_t ofsj5 = (pcgpi->ofs * sizeof(UINT32)) / sizeof(Jvalue5);
    size_t celj5 = (pcgpi->cel * sizeof(UINT32)) / sizeof(Jvalue5);

    // start at the end of the current partition
    Jvalue5* pJ0 = reinterpret_cast<Jvalue5*>(_pHostBuf) + ofsj5 + celj5 - 1;
    Jvalue5* pJlimit = pJ0 + min2(maxnJ, pcgpi[1].cel);

    // start with a value in the J list at the end of the partition
    UINT64 vPrev = *reinterpret_cast<UINT64*>(pJ0) & j5mask;

    // scan forward to find the first out-of-order J value in the table
    Jvalue5* pJ = pJ0;
    while( ++pJ <= pJlimit )
    {
        UINT64 v = *reinterpret_cast<UINT64*>(pJ) & j5mask;
        if( v < vPrev )
        {
            /* At this point we have discovered the first 5-byte value in the table that is not in the initial J list. */
            break;
        }

        /* At this point we are either still in the same J list or in a subsequent J list whose first J value happens
            to exceed that of the last J value in the previous list.  Either way, we want to keep iterating. */
        vPrev = v;
    }

    /* At this point pJ points to the first out-of-order J value. */

    // save the new J table partition size (including one extra UINT32 of padding)
    celj5 = pJ - reinterpret_cast<Jvalue5*>(_pHostBuf);         // number of Jvalue5 values in the table
    size_t celui4 = (celj5 * sizeof(Jvalue5)) / sizeof(UINT32); // number of UINT32 values in the table
    return max2(pcgpi->cel, celui4-pcgpi->ofs) + 1;             // number of UINT32 values in the partition data table
}

/// [private] method loadLUT
void tuGpu::loadLUT( GpuInfo* _pgi, CudaGmemLayout* _pcgl )
{
    CRVALIDATOR;

    try
    {
        UINT32* pHostBuf = reinterpret_cast<UINT32*>(_pcgl->HostBuf->GetBufferPointer());
        CudaGmemLayout::CudaGmemPartitionInfo* pcgpi = _pcgl->part + m_gpuDeviceOrdinal;

        if( _pcgl->layoutType == cgltPinned )
        {
            // use the table in page-locked host memory (CUDA-addressable "pinned" memory)
            pcgpi->p = pHostBuf;

            CDPrint( cdpCD2, "%s: initialized %s buffer (%llu bytes) at 0x%016llx for device %d (ordinal %d) in CUDA pinned memory", __FUNCTION__,
                                _pcgl->tag, pcgpi->cel*sizeof(UINT32), pcgpi->p, _pgi->deviceId, m_gpuDeviceOrdinal );
        }
        else
        {
            /* Initialize the CUDA global memory buffer on the current GPU device (see AriocP::loadHJ).

               If the buffer is partitioned across multiple GPUs, each buffer's contents must
                overlap the subsequent buffer by one UINT32 (see tuSetupN10::getV5() for example).
                We deal with this by allocating and copying one additional UINT32, and we clamp
                that allocation so that it does not extend past the end of the data.
            */
            size_t cel = pcgpi->cel;

            if( _pcgl->layoutType == cgltGPmem )
            {
                // enable "peer access" for the current GPU device
                _pgi->pCDB->SetPeerMemoryAccess( true );

                if( m_gpuDeviceOrdinal < (_pcgl->nPartitions-1) )
                {
                    /* Ensure that each partition overlaps the subsequent partition:
                        - for an H table: one 32-bit value is sufficient because each H value is 5 bytes in size
                        - for a J table: we sniff out the start of the last J list in the partition and ensure that
                           all of its J values are present
                    */
                    switch( _pcgl->tag[0] )
                    {
                        case 'H':
                            ++cel;
                            break;

                        case 'J':
                            cel = sniffJeol( _pcgl, pHostBuf );
                            CDPrint( cdpCDb, "%s: sniffJeol for %s returns %lld (partition cel=%lld)", __FUNCTION__, _pcgl->tag, cel, pcgpi->cel );
                            break;

                        default:
                            throw new ApplicationException( __FILE__, __LINE__, "unexpected GPU layout tag '%s'", _pcgl->tag );
                    }
                }
            }

            // allocate a CUDA global memory buffer to contain the table partition
            if( pcgpi->PartBuf )
                throw new ApplicationException( __FILE__, __LINE__, "%s: cannot reallocate CUDA global memory buffer for %s", __FUNCTION__, _pcgl->tag );
            CREXEC( pcgpi->PartBuf = new CudaGlobalPtr<UINT32>( cel, true ) );

            // copy the table partition data from the host buffer to the CUDA global memory buffer
            CREXEC( pcgpi->PartBuf->CopyToDevice( pHostBuf+pcgpi->ofs, cel ) );
            pcgpi->p = pcgpi->PartBuf->p;

            CDPrint( cdpCD2, "%s: initialized %s buffer (%llu bytes) at 0x%016llx for device %d (ordinal %d) in CUDA global memory", __FUNCTION__,
                                _pcgl->tag, pcgpi->PartBuf->cb, pcgpi->p, _pgi->deviceId, m_gpuDeviceOrdinal );

            // discard the host buffer if it is no longer needed
            if( InterlockedIncrement( &_pcgl->HostBuf->n ) == static_cast<UINT32>(m_pab->nGPUs) )
            {
                _pcgl->HostBuf->Free();
                CDPrint( cdpCD2, "%s: freed %s staging buffer", __FUNCTION__, _pcgl->tag );
                _pcgl->HostBuf = NULL;
            }
        }
    }
    catch( ApplicationException* _pex )
    {
        _pex->SetCallerExceptionInfo( __FILE__, __LINE__, "unable to allocate CUDA-addressable memory for %s table", _pcgl->tag );
        throw _pex;
    }
}

/// [private] method unloadLUT
void tuGpu::unloadLUT( CudaGmemLayout* _pcgl )
{
    CudaGmemLayout::CudaGmemPartitionInfo* pcgpi = _pcgl->part + m_gpuDeviceOrdinal;
    if( pcgpi->PartBuf )
    {
        pcgpi->PartBuf->Free();
        delete pcgpi->PartBuf;
        pcgpi->PartBuf = NULL;
    }
}
#pragma endregion

#pragma region protected methods
/// [protected] method loadR
void tuGpu::loadR( GpuInfo* pgi )
{
    if( m_pab->usePinnedR )
        pgi->pR = m_pab->R;
    else
    {
        // copy R table to CUDA global memory
        CRVALIDATOR;
        SET_CUDAGLOBALPTR_TAG( pgi->bufR, "R" );
        CREXEC( pgi->bufR.Realloc( m_pab->celR, false ) );
        CREXEC( pgi->bufR.CopyToDevice( m_pab->R, pgi->bufR.Count ) );
        pgi->pR = pgi->bufR.p;

        CDPrint( cdpCD2, "%s: initialized R buffer (%llu bytes) at 0x%016llx in CUDA global memory", __FUNCTION__, pgi->bufR.cb, pgi->pR );
    }
}

/// [protected] method loadHJ
void tuGpu::loadHJ( GpuInfo* pgi )
{
    // load lookup tables
    loadLUT( pgi, &m_pab->cglHn );
    loadLUT( pgi, &m_pab->cglJn );
    loadLUT( pgi, &m_pab->cglHg );
    loadLUT( pgi, &m_pab->cglJg );
}

/// [protected] method loadQ
void tuGpu::loadQ( QBatch* pqb )
{
    CRVALIDATOR;

    HiResTimer hrt(us);

    try
    {
        /* CUDA global memory layout for all kernels:

             low:   Qw      Qwarps
                    Qi      interleaved Q sequence data
                    ...     ...

             high:  ...     ...
        */

#if TODO_CHOP_WHEN_DEBUGGED
        pqb->pgi->pCGA->DumpUsage( __FUNCTION__ );
#endif

        /* Load the Qwarp and Q sequence data into global memory on the current GPU; these allocations are freed just before the
            current QBatch instance is released (see tuTailP::main) */

        // Qwarp data
        CREXEC( pqb->DB.Qw.Alloc( cgaLow, pqb->QwBuffer.n, false ) );
        CREXEC( pqb->DB.Qw.CopyToDevice( pqb->QwBuffer.p, pqb->DB.Qw.Count ) );
        pqb->DB.Qw.n = pqb->QwBuffer.n;
        SET_CUDAGLOBALPTR_TAG( pqb->DB.Qw, "DB.Qw" );

        // interleaved Q sequence data
        CREXEC( pqb->DB.Qi.Alloc( cgaLow, pqb->QiBuffer.n, false ) );
        CREXEC( pqb->DB.Qi.CopyToDevice( pqb->QiBuffer.p, pqb->DB.Qi.Count ) );
        pqb->DB.Qi.n = pqb->QiBuffer.n;
        SET_CUDAGLOBALPTR_TAG( pqb->DB.Qi, "DB.Qi" );
    }
    catch( ApplicationException* pex )
    {
        CRTHROW;
    }

    // performance metrics
    InterlockedExchangeAdd( &AriocBase::aam.us.XferQ, hrt.GetElapsed(false) );
}

/// [protected] method unloadQ
void tuGpu::unloadQ( QBatch* pqb )
{
    // free the GPU buffers that contain the Qwarp and interleaved Q sequence data   
    pqb->DB.Qi.Free();
    pqb->DB.Qw.Free();

#if TODO_CHOP_WHEN_DEBUGGED
    CDPrint( cdpCD0, "tuGpu::unloadQ dumps the global allocation for batch 0x%016llx...", pqb );
    pqb->pgi->pCGA->DumpUsage();
#endif
}

/// [protected] method unloadHJ
void tuGpu::unloadHJ( GpuInfo* _pgi )
{
    HiResTimer hrt( ms );

    // delete CUDA global memory buffers
    unloadLUT( &m_pab->cglHn );
    unloadLUT( &m_pab->cglJn );
    unloadLUT( &m_pab->cglHg );
    unloadLUT( &m_pab->cglJg );

    // conditionally disable "peer access" for the current GPU device
    if( (m_pab->cglHn.layoutType == cgltGPmem) || (m_pab->cglJn.layoutType == cgltGPmem) || (m_pab->cglHg.layoutType == cgltGPmem) || (m_pab->cglJg.layoutType == cgltGPmem) )
        _pgi->pCDB->SetPeerMemoryAccess( false );

    // performance metrics
    AriocBase::aam.ms.UnloadHJ += hrt.GetElapsed( false );
}
#pragma endregion
