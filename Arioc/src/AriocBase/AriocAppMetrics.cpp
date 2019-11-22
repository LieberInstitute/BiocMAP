/*
   AriocAppMetrics.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#include <thrust/version.h>

#pragma region static member variable definitions
const char* AriocAppMetrics::m_fmt_n1 = "%s%-*s: %*lld";
const char* AriocAppMetrics::m_fmt_n2 = "%s%-*s: %*lld (%lld SAM records)";
const char* AriocAppMetrics::m_fmt_n3 = "%s%-*s: %*u";
const char* AriocAppMetrics::m_fmt_cb1 = "%s%-*s: %*lld (%2.1fGB)";
const char* AriocAppMetrics::m_fmt_cb2 = "%s%-*s: %*lld (%2.1fGB) (max %u J/seed, max subId %d)";
const char* AriocAppMetrics::m_fmt_cb3 = "%s%-*s: %*lld (%2.1fGB) %s";
const char* AriocAppMetrics::m_fmt_ms1 = "%s%-*s: %*lld";
#pragma endregion

#pragma region AriocAppMetrics: constructor/destructor
/// default constructor
AriocAppMetrics::AriocAppMetrics(): GPUmask(0)
{
}

/// destructor
AriocAppMetrics::~AriocAppMetrics()
{
}
#pragma endregion

#pragma region AriocAppMetrics: protected methods
/// [protected] method emitGpuInfo
void AriocAppMetrics::emitGpuInfo( CDPrintFlags cdp, const char* indent2, INT32 cchSpan1 )
{
    const char* fmtGmemG = "%lld (%2.1fGB)";
    const char* fmtGmemP = "%lld (%2.1f%%)";

    // emit the CUDA runtime version
    CDPrint( cdp, "%s%-*s: %s", indent2, cchSpan1, "CUDA version", CudaCommon::GetCudaVersionString() );

    // emit the CUDA Thrust library version
    CDPrint( cdp, "%s%-*s: %d.%d.%d", indent2, cchSpan1, "Thrust library version", THRUST_MAJOR_VERSION, THRUST_MINOR_VERSION, THRUST_SUBMINOR_VERSION );

    /* GPU device info is emitted in a format like this:

           GPU resources                  : global memory         ECC  threads  CGA total             CGA used
            0: Tesla V100-SXM2-32GB       : 34058272768 (31.7GB)  on   163840   12352290816 (12.1GB)  12352290816 (100.0%)
            1: TITAN Xp                   : (not used)
            2: Tesla M40                  :  4283142144 ( 4.1GB)   4283142144 ( 0.1%)
    */

    this->n.GPUs = CudaCommon::GetDeviceCount();

    // compute column widths
    size_t cchGmem = 0;
    size_t cchCGAtotal = 0;
    size_t cchCGAused = 0;
    for( INT16 deviceId=0; deviceId<this->n.GPUs; ++deviceId )
    {
        cudaDevicePropEx cdpx;
        CudaCommon::GetDeviceProperties( &cdpx, deviceId );

        char s[32];
        sprintf_s( s, sizeof s, fmtGmemG, cdpx.totalGlobalMem, cdpx.totalGlobalMem/m_bytesPerGB );
        cchGmem = max2(cchGmem, strlen( s ));

        sprintf_s( s, sizeof s, fmtGmemG, this->cb.TotalCGA[deviceId], this->cb.TotalCGA[deviceId]/m_bytesPerGB );
        cchCGAtotal = max2( cchCGAtotal, strlen( s ) );

        sprintf_s( s, sizeof s, fmtGmemP, this->cb.UsedCGA[deviceId], 100.0*this->cb.UsedCGA[deviceId]/this->cb.TotalCGA[deviceId] );
        cchCGAused = max2( cchCGAused, strlen( s ) );
    }

    // display the table
    CDPrint( cdp, "%s%-*s: %-*s  %s  %-*s  %-*s", indent2, cchSpan1, "GPU resources", cchGmem, "global memory", "ECC  threads", cchCGAtotal, "CGA total", cchCGAused, "CGA used" );
    for( INT16 deviceId=0; deviceId<this->n.GPUs; ++deviceId )
    {
        cudaDevicePropEx cdpx;
        CudaCommon::GetDeviceProperties( &cdpx, deviceId );

        char banner[64];
        sprintf_s( banner, sizeof banner, " %2d: %s", deviceId, cdpx.name );

        if( this->GPUmask & (1 << deviceId) )
        {
            char sGmem[32];
            sprintf_s( sGmem, sizeof sGmem, fmtGmemG, cdpx.totalGlobalMem, cdpx.totalGlobalMem/m_bytesPerGB );
            const char* ecc = cdpx.ECCEnabled ? "on" : "off";

            INT32 threads = cdpx.maxThreadsPerMultiProcessor * cdpx.multiProcessorCount;

            char sCGAtotal[32];
            sprintf_s( sCGAtotal, sizeof sCGAtotal, fmtGmemG, this->cb.TotalCGA[deviceId], this->cb.TotalCGA[deviceId]/m_bytesPerGB );

            char sCGAused[32];
            sprintf_s( sCGAused, sizeof sCGAused, fmtGmemP, this->cb.UsedCGA[deviceId], 100.0*this->cb.UsedCGA[deviceId]/this->cb.TotalCGA[deviceId] );

            CDPrint( cdp, "%s%-*s: %-*s  %-3s  %-7d  %-*s  %-*s", indent2, cchSpan1, banner, cchGmem, sGmem, ecc, threads, cchCGAtotal, sCGAtotal, cchCGAused, sCGAused );
        }
        else
            CDPrint( cdp, "%s%-*s: (not used)", indent2, cchSpan1, banner );
    }
}

/// [protected] method emitGpuMemInfo
void AriocAppMetrics::emitGpuMemInfo( CDPrintFlags cdp, const char* indent2, INT32 cchSpan1, INT32 cchSpan2 )
{
    // compute total GPU memory used for R and HJ data
    UINT64 cbRHJgmem = this->cb.Rgmem + (this->cb.HJpartitioned ? 0 : this->cb.HJgmem);
    UINT64 cbHJpartitioned = this->cb.HJpartitioned ? this->cb.HJgmem : 0;

    // compute the column width for numerical data
    if( cchSpan2 == -1 )
    {
        // find the maximum numerical value
        UINT64 maxval = 0;
        if( this->cb.Rpinned > maxval )         maxval = this->cb.Rpinned;
        if( this->cb.HJpinned > maxval )        maxval = this->cb.HJpinned;
        if( cbRHJgmem > maxval )                maxval = cbRHJgmem;
        if( cbHJpartitioned > maxval )          maxval = cbHJpartitioned;

        cchSpan2 = maxval ? static_cast<INT32>(log10( static_cast<double>(maxval) )) + 1 : 1;
    }

    CDPrint( cdp, "%sLUT memory resources:", indent2 );
    if( this->cb.Rpinned || this->cb.HJpinned )
    {
        CDPrint( cdp, m_fmt_cb1, indent2, cchSpan1, " pinned (page-locked) RAM", cchSpan2, this->cb.HJpinned, this->cb.HJpinned/m_bytesPerGB );
        if( this->cb.Rpinned )
            CDPrint( cdp, m_fmt_cb1, indent2, cchSpan1, "  R sequence data", cchSpan2, this->cb.Rpinned, this->cb.Rpinned/m_bytesPerGB );
        if( this->cb.Hnp )
            CDPrint( cdp, m_fmt_cb1, indent2, cchSpan1, "  H lookup table (nongapped)", cchSpan2, this->cb.Hnp, this->cb.Hnp/m_bytesPerGB );
        if( this->cb.Hgp )
            CDPrint( cdp, m_fmt_cb1, indent2, cchSpan1, "  H lookup table (gapped)", cchSpan2, this->cb.Hgp, this->cb.Hgp/m_bytesPerGB );
        if( this->cb.Jnp )
            CDPrint( cdp, m_fmt_cb2, indent2, cchSpan1, "  J lookup table (nongapped)", cchSpan2, this->cb.Jnp, this->cb.Jnp/m_bytesPerGB, this->n.maxnJn, this->n.maxSubId );
        if( this->cb.Jgp )
            CDPrint( cdp, m_fmt_cb2, indent2, cchSpan1, "  J lookup table (gapped)", cchSpan2, this->cb.Jgp, this->cb.Jgp/m_bytesPerGB, this->n.maxnJg, this->n.maxSubId );
    }
    if( cbRHJgmem )
    {
        CDPrint( cdp, m_fmt_cb1, indent2, cchSpan1, " CUDA global memory (per GPU)", cchSpan2, cbRHJgmem, cbRHJgmem/m_bytesPerGB );

        if( this->cb.Rgmem )
            CDPrint( cdp, m_fmt_cb1, indent2, cchSpan1, "  R sequence data", cchSpan2, this->cb.Rgmem, this->cb.Rgmem/m_bytesPerGB );

        if( cbHJpartitioned )
            CDPrint( cdp, m_fmt_cb1, indent2, cchSpan1, " CUDA global memory (P2P)", cchSpan2, cbHJpartitioned, cbHJpartitioned/m_bytesPerGB );

        if( this->cb.Hng )
            CDPrint( cdp, m_fmt_cb1, indent2, cchSpan1, "  H lookup table (nongapped)", cchSpan2, this->cb.Hng, this->cb.Hng/m_bytesPerGB );
        if( this->cb.Hgg )
            CDPrint( cdp, m_fmt_cb1, indent2, cchSpan1, "  H lookup table (gapped)", cchSpan2, this->cb.Hgg, this->cb.Hgg/m_bytesPerGB );
        if( this->cb.Jng )
            CDPrint( cdp, m_fmt_cb2, indent2, cchSpan1, "  J lookup table (nongapped)", cchSpan2, this->cb.Jng, this->cb.Jng/m_bytesPerGB, this->n.maxnJn, this->n.maxSubId );
        if( this->cb.Jgg )
            CDPrint( cdp, m_fmt_cb2, indent2, cchSpan1, "  J lookup table (gapped)", cchSpan2, this->cb.Jgg, this->cb.Jgg/m_bytesPerGB, this->n.maxnJg, this->n.maxSubId );
    }
}

/// [protected] method emitCpuInfo
void AriocAppMetrics::emitCpuInfo( CDPrintFlags cdp, const char* indent2, INT32 cchSpan1, INT32 cchSpan2 )
{
    CpuInfo cpui;

    // compute the column width for numerical data
    if( cchSpan2 == -1 )
    {
        // find the maximum numerical value
        UINT64 maxval = 0;
        if( this->cb.AvailableCPU > maxval )    maxval = this->cb.AvailableCPU;

        cchSpan2 = maxval ? static_cast<INT32>(log10( static_cast<double>(maxval) )) + 1 : 1;
    }

    CDPrint( cdp, "%sHost (CPU) resources:", indent2 );
    if( cpui.Description[0] )
        CDPrint( cdp, "%s%-*s: %s", indent2, cchSpan1, " CPU description", cpui.Description );
    CDPrint( cdp, "%s%-*s: %d", indent2, cchSpan1, " concurrent threads", this->n.CPUThreads );
    CDPrint( cdp, m_fmt_cb1, indent2, cchSpan1, " available physical RAM", cchSpan2, this->cb.AvailableCPU, this->cb.AvailableCPU/m_bytesPerGB );
}

/// [protected] method emitSamInfoU
void AriocAppMetrics::emitSamInfoU( CDPrintFlags cdp, const char* indent2, INT32 cchSpan1, INT32 cchSpan2 )
{
    CDPrint( cdp, "%sData transfers:", indent2 );
    AriocBase::aam.PrintPipelineInfo( cdp, 1, static_cast<INT32>(strlen( indent2 )), cchSpan1 - 1, cchSpan2 );

    // compute some additional metrics
    UINT64 nTotalMapped = this->n.u.nReads_Mapped1 + this->n.u.nReads_Mapped2;
    UINT64 nTotalSAMmapped = this->sam.u.nReads_Mapped1 + this->sam.u.nReads_Mapped2;
    UINT64 nTotalSAM = nTotalSAMmapped + this->sam.u.nReads_Unmapped;
    this->n.pctReadsMapped = 100.0 * nTotalMapped / this->sam.nCandidateQ;

    // compute the column width for numerical data
    if( cchSpan2 == -1 )
    {
        // find the maximum numerical value
        UINT64 maxval = this->sam.nCandidateQ;
        if( this->n.u.nReads_Mapped1 > maxval )     maxval = this->n.u.nReads_Mapped1;
        if( this->n.u.nReads_Mapped2 > maxval )     maxval = this->n.u.nReads_Mapped2;
        if( this->sam.u.nReads_Unmapped > maxval)   maxval = this->sam.u.nReads_Unmapped;
        if( this->n.Nmax > maxval)                  maxval = this->n.Nmax;
        if( this->n.bwmax > maxval)                 maxval = this->n.bwmax;
        if( this->n.SketchBits > maxval )           maxval = this->n.SketchBits;

        cchSpan2 = maxval ? static_cast<INT32>(log10( static_cast<double>(maxval) )) + 1 : 1;
    }

    CDPrint( cdp, "%sSAM output:", indent2 );
    CDPrint( cdp, m_fmt_n2, indent2, cchSpan1, " reads", cchSpan2, this->sam.nCandidateQ, nTotalSAM );
    CDPrint( cdp, "%s%-*s: %*lld (%4.2f%%)", indent2, cchSpan1, " mapped reads", cchSpan2, nTotalMapped, this->n.pctReadsMapped );
    CDPrint( cdp, m_fmt_n2, indent2, cchSpan1, "  with 1 mapping",          cchSpan2, this->n.u.nReads_Mapped1, this->sam.u.nReads_Mapped1 );
    CDPrint( cdp, m_fmt_n2, indent2, cchSpan1, "  with 2 or more mappings", cchSpan2, this->n.u.nReads_Mapped2, this->sam.u.nReads_Mapped2 );
    UINT64 nTotalUnmapped = this->sam.nCandidateQ - nTotalMapped;
    CDPrint( cdp, m_fmt_n2, indent2, cchSpan1, " unmapped reads",           cchSpan2, nTotalUnmapped, this->sam.u.nReads_Unmapped );

    if( this->n.SketchBits )
    {
        CDPrint( cdp, "%sKMH output:", indent2 );
        CDPrint( cdp, m_fmt_n3, indent2, cchSpan1, " sketch bits", cchSpan2, this->n.SketchBits );
    }
}

/// [protected] method emitSamInfoP
void AriocAppMetrics::emitSamInfoP( CDPrintFlags cdp, const char* indent2, INT32 cchSpan1, INT32 cchSpan2 )
{
    CDPrint( cdp, "%sData transfers:", indent2 );
    AriocBase::aam.PrintPipelineInfo( cdp, 1, static_cast<INT32>(strlen( indent2 )), cchSpan1 - 1, cchSpan2 );

    // compute some additional metrics
    UINT64 nCandidatePairs = this->sam.nCandidateQ / 2;
    UINT64 nPairsConcordant = this->sam.p.nPairs_Concordant1 + this->sam.p.nPairs_Concordant2;
    UINT64 nMatesInUnmappedPairs = this->sam.p.nMates_Unmapped + this->sam.p.nMates_MappedU1 + this->sam.p.nMates_MappedU2;
    UINT64 nTotalMappedMates = 2 * (nPairsConcordant + this->sam.p.nPairs_Discordant) + this->sam.p.nMates_MappedU1 + this->sam.p.nMates_MappedU2;
    this->n.pctReadsMapped = 100.0 * nTotalMappedMates / this->sam.nCandidateQ;
    UINT64 nTotal = this->sam.p.nRowsConcordant + this->sam.p.nRowsDiscordant + this->sam.p.nRowsRejected + this->sam.p.nRowsUnmapped;
    double pctConcordant = 100.0 * nPairsConcordant / nCandidatePairs;

    // compute the column width for numerical data
    if( cchSpan2 == -1 )
    {
        // find the maximum numerical value
        UINT64 maxval = 0;
        if( nCandidatePairs > maxval )              maxval = nCandidatePairs;
        if( this->sam.p.nRowsConcordant > maxval )  maxval = this->sam.p.nRowsConcordant;
        if( this->sam.p.nRowsDiscordant > maxval)   maxval = this->sam.p.nRowsDiscordant;
        if( this->sam.p.nRowsRejected > maxval)     maxval = this->sam.p.nRowsRejected;
        if( this->sam.p.nRowsUnmapped > maxval)     maxval = this->sam.p.nRowsUnmapped;
        if( nMatesInUnmappedPairs > maxval)         maxval = nMatesInUnmappedPairs;
        if( nTotalMappedMates > maxval)             maxval = nTotalMappedMates;
        if( this->n.DuplicateMappings > maxval )    maxval = this->n.DuplicateMappings;
        if( this->n.Nmax > maxval )                 maxval = this->n.Nmax;
        if( this->n.bwmax > maxval)                 maxval = this->n.bwmax;
        if( this->n.SketchBits > maxval)            maxval = this->n.SketchBits;

        cchSpan2 = maxval ? static_cast<INT32>(log10( static_cast<double>(maxval) )) + 1 : 1;
    }

    CDPrint( cdp, "%sSAM output:", indent2 );
    CDPrint( cdp, m_fmt_n2, indent2, cchSpan1, " pairs", cchSpan2, nCandidatePairs, nTotal );
            
    CDPrint( cdp, "%s%-*s: %*lld (%lld SAM records) (%4.2f%%)", indent2, cchSpan1, "  concordant pairs", cchSpan2, nPairsConcordant, this->sam.p.nRowsConcordant, pctConcordant );
    CDPrint( cdp, m_fmt_n1, indent2, cchSpan1, "   with 1 mapping", cchSpan2, this->sam.p.nPairs_Concordant1 );
    CDPrint( cdp, m_fmt_n1, indent2, cchSpan1, "   with 2 or more mappings", cchSpan2, this->sam.p.nPairs_Concordant2 );

    CDPrint( cdp, m_fmt_n2, indent2, cchSpan1, "  discordant pairs", cchSpan2, this->sam.p.nPairs_Discordant, this->sam.p.nRowsDiscordant );
    CDPrint( cdp, m_fmt_n2, indent2, cchSpan1, "  rejected pairs", cchSpan2, this->sam.p.nPairs_Rejected, this->sam.p.nRowsRejected );
    CDPrint( cdp, m_fmt_n2, indent2, cchSpan1, "  unmapped pairs", cchSpan2, this->sam.p.nPairs_Unmapped, this->sam.p.nRowsUnmapped );
    
    CDPrint( cdp, m_fmt_n1, indent2, cchSpan1, "  mates not in paired mappings", cchSpan2, nMatesInUnmappedPairs );
    CDPrint( cdp, m_fmt_n1, indent2, cchSpan1, "   with no mappings", cchSpan2, this->sam.p.nMates_Unmapped );
    CDPrint( cdp, m_fmt_n1, indent2, cchSpan1, "   with 1 mapping", cchSpan2, this->sam.p.nMates_MappedU1 );
    CDPrint( cdp, m_fmt_n1, indent2, cchSpan1, "   with 2 or more mappings", cchSpan2, this->sam.p.nMates_MappedU2 );

    CDPrint( cdp, "%s%-*s: %*lld (%4.2f%%)", indent2, cchSpan1, " total mapped mates", cchSpan2, nTotalMappedMates, this->n.pctReadsMapped );
    CDPrint( cdp, m_fmt_n1, indent2, cchSpan1, " duplicate mappings", cchSpan2, this->n.DuplicateMappings );
    CDPrint( cdp, m_fmt_n3, indent2, cchSpan1, " maximum Q length", cchSpan2, this->n.Nmax );
    CDPrint( cdp, m_fmt_n3, indent2, cchSpan1, " maximum diagonal band width", cchSpan2, this->n.bwmax );
    CDPrint( cdp, "%s%-*s: %*.1f (%2.1f)", indent2, cchSpan1, " TLEN mean (stdev)", cchSpan2, this->n.p.TLEN_mean, this->n.p.TLEN_sd );

    if( this->n.SketchBits )
    {
        CDPrint( cdp, "%sKMH output:", indent2 );
        CDPrint( cdp, m_fmt_n3, indent2, cchSpan1, " sketch bits", cchSpan2, this->n.SketchBits );
    }
}

/// [protected] method emitSamInfo
void AriocAppMetrics::emitSamInfo( CDPrintFlags cdp, const char* indent2, INT32 cchSpan1, INT32 cchSpan2 )
{
    if( this->sam.p.nPairs_Concordant1 || this->sam.p.nPairs_Discordant || this->sam.p.nPairs_Rejected || this->sam.p.nPairs_Unmapped )
        emitSamInfoP( cdp, indent2, cchSpan1, cchSpan2 );
    else
    if( this->sam.u.nReads_Mapped1 || this->sam.u.nReads_Mapped2 || this->sam.u.nReads_Unmapped )
        emitSamInfoU( cdp, indent2, cchSpan1, cchSpan2 );
}

/// [protected] method emitSummaryInfo
void AriocAppMetrics::emitSummaryInfo( CDPrintFlags cdp, const char* indent2, INT32 cchSpan1, INT32 cchSpan2 )
{
    // compute some additional metrics
    UINT64 msLUTs = this->ms.LoadR + this->ms.LoadHJ + this->ms.UnloadR + this->ms.UnloadHJ;
    UINT64 msOverhead = this->ms.InitializeGPUs + this->ms.SniffQ + msLUTs;
    if( this->ms.App > msOverhead )
        this->ms.Aligners = this->ms.App - msOverhead;
    else
        this->ms.Aligners = 1;  // (can't determine elapsed time for aligners at millisecond resolution, so report 1ms instead)

    // compute the column width for numerical data
    if( cchSpan2 == -1 )
    {
        // find the maximum numerical value
        UINT64 maxval = this->ms.App;

        cchSpan2 = maxval ? static_cast<INT32>(log10( static_cast<double>(maxval) )) + 1 : 1;
    }

    CDPrint( cdp, "%sBatches:", indent2 );
    CDPrint( cdp, m_fmt_n3, indent2, cchSpan1, " instances", cchSpan2, this->n.BatchInstances );
    CDPrint( cdp, m_fmt_n3, indent2, cchSpan1, " completed", cchSpan2, this->n.BatchesCompleted );
    CDPrint( cdp, "%sElapsed (ms):", indent2 );
    CDPrint( cdp, m_fmt_ms1, indent2, cchSpan1, " total", cchSpan2, this->ms.App );
    CDPrint( cdp, m_fmt_ms1, indent2, cchSpan1, " initialize GPUs", cchSpan2, this->ms.InitializeGPUs );
    CDPrint( cdp, m_fmt_ms1, indent2, cchSpan1, " partition Q files", cchSpan2, this->ms.SniffQ );
    CDPrint( cdp, m_fmt_ms1, indent2, cchSpan1, " load/unload R,H,J", cchSpan2, msLUTs );
    CDPrint( cdp, "%s%-*s: %*lld (%1.0f Q/second)", indent2, cchSpan1, " aligners", cchSpan2, this->ms.Aligners, m_msPerSec*this->n.Q/this->ms.Aligners );
}
#pragma endregion

#pragma region AriocAppMetrics: public methods
/// <summary>
/// Displays or prints the non-null contents of the current <c>AriocTaskUnitMetrics</c> instance.
/// </summary>
/// <param name="cdp">flags for <c>CDPrint</c></param>
/// <param name="flags">specifies which metric(s) are emitted even when zero</param>
/// <param name="description">description header</param>
/// <param name="cchIndent1">level 1 indentation (for description header only)</param>
/// <param name="cchIndent2">level 2 indentation (for labels; includes level 1 indentation)</param>
/// <param name="cchSpan1">column width for labels</param>
/// <param name="cchSpan2">[optional] column width for right-aligned numerical values</param>
void AriocAppMetrics::Print( CDPrintFlags cdp, const char* description, INT32 cchIndent1, INT32 cchIndent2, INT32 cchSpan1, INT32 cchSpan2 )
{
    // set up two strings of spaces to serve as indentation
    WinGlobalPtr<char> indent1(cchIndent1+1, false);
    memset( indent1.p, ' ', cchIndent1 );
    indent1.p[cchIndent1] = 0;

    cchIndent2 += cchIndent1;       // (second-level indentation is relative to the first-level indentation)
    WinGlobalPtr<char> indent2(cchIndent2+1, false);
    memset( indent2.p, ' ', cchIndent2 );
    indent2.p[cchIndent2] = 0;

    // use the description (if any) as a header
    if( description )
        CDPrint( cdp, "%s%s:", indent1.p, description );

    /* resources */
    emitGpuInfo( cdp, indent2.p, cchSpan1 );
    emitCpuInfo( cdp, indent2.p, cchSpan1, cchSpan2 );
    emitGpuMemInfo( cdp, indent2.p, cchSpan1, cchSpan2 );

    /* mappings */
    emitSamInfo( cdp, indent2.p, cchSpan1, cchSpan2 );

    /* summary */
    emitSummaryInfo( cdp, indent2.p, cchSpan1, cchSpan2 );
}

/// <summary>
/// Displays or prints the non-null contents of the &quot;pipeline info&quot; members of the <c>AriocAppMetrics</c> struct
/// </summary>
void AriocAppMetrics::PrintPipelineInfo( CDPrintFlags cdp, INT32 cchIndent1, INT32 cchIndent2, INT32 cchSpan1, INT32 cchSpan2 )
{
    // set up two strings of spaces to serve as indentation
    WinGlobalPtr<char> indent1(cchIndent1+1, false);
    memset( indent1.p, ' ', cchIndent1 );
    indent1.p[cchIndent1] = 0;

    cchIndent2 += cchIndent1;       // (second-level indentation is relative to the first-level indentation)
    WinGlobalPtr<char> indent2(cchIndent2+1, false);
    memset( indent2.p, ' ', cchIndent2 );
    indent2.p[cchIndent2] = 0;

    // convert microseconds to milliseconds
    UINT64 msXferQ = usToMs( this->us.XferQ );
    UINT64 msXferR = usToMs( this->us.XferR );
    UINT64 msXferJ = usToMs( this->us.XferJ );
    UINT64 msProbeHJ = usToMs( this->us.ProbeHJ );

    // compute the column width for numerical data
    if( cchSpan2 == -1 )
    {
        // find the maximum numerical value
        UINT64 maxval = this->ms.LoadQ;
        if( this->ms.BuildQwarps > maxval )     maxval = this->ms.BuildQwarps;
        if( msXferQ > maxval )                  maxval = msXferQ;
        if( msXferR > maxval )                  maxval = msXferR;
        if( msXferJ > maxval )                  maxval = msXferJ;
        if( this->ms.XferMappings > maxval )    maxval = this->ms.XferMappings;
        if( msProbeHJ > maxval )                maxval = msProbeHJ;
        if( this->ms.WriteA > maxval )          maxval = this->ms.WriteA;

        // compute the number of digits in the maximum numerical value
        cchSpan2 = maxval ? static_cast<INT32>(log10( static_cast<double>(maxval) )) + 1 : 1;
    }

    CDPrint( cdp, m_fmt_ms1, indent2.p, cchSpan1, "disk --> CPU: load Q data", cchSpan2, this->ms.LoadQ );
    CDPrint( cdp, m_fmt_ms1, indent2.p, cchSpan1, "CPU: build Qwarps", cchSpan2, this->ms.BuildQwarps );
    CDPrint( cdp, m_fmt_ms1, indent2.p, cchSpan1, "CPU --> GPU: transfer Q data", cchSpan2, msXferQ );
    CDPrint( cdp, m_fmt_ms1, indent2.p, cchSpan1, "CPU --> GPU: transfer R data", cchSpan2, msXferR );
    CDPrint( cdp, m_fmt_ms1, indent2.p, cchSpan1, "GPU --> CPU: transfer J data", cchSpan2, msXferJ );
    CDPrint( cdp, m_fmt_ms1, indent2.p, cchSpan1, "GPU --> CPU: transfer mappings", cchSpan2, this->ms.XferMappings );
    CDPrint( cdp, m_fmt_ms1, indent2.p, cchSpan1, "GPU --> GPU: probe H and J LUTs", cchSpan2, msProbeHJ );
    CDPrint( cdp, m_fmt_ms1, indent2.p, cchSpan1, "CPU --> disk: write alignments", cchSpan2, this->ms.WriteA );
}
#pragma endregion
