/*
  AriocBase.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#pragma once
#define __AriocBase__

#pragma region enums
enum BRLEAbyteType : UINT8
{
    bbMatch     =   0x00,   // 00
    bbGapQ      =   0x01,   // 01 (deletion from R)
    bbMismatch  =   0x02,   // 10
    bbGapR      =   0x03    // 11 (insertion into R)
};
#pragma endregion

#pragma region structs
#pragma pack(push, 1)
/* A BRLEA (binary run-length encoded alignment) is a string of bytes packed as follows:

                        7  6  5  4  3  2  1  0
    match               0  0  |--- cchRun ---|
    deletion (from R)   0  1  |--- cchRun ---|
    mismatch            1  0  |--- cchRun ---|
    insert (Q into R)   1  1  |--- cchRun ---|
*/
struct BRLEAbyte
{
    UINT8   bbRunLength : 6;    // 0..5: run length
    UINT8   bbType      : 2;    // 6..7: BRLEAbyteType
};
#define packBRLEAbyte(T,R)          (static_cast<UINT8>((T << 6) | (R)))

struct BRLEAheader
{
    UINT32  qid;    // packed reference to the Qwarp for the Q sequence
    INT16   V;      // V score
    INT16   Ma;     // number of R symbols spanned by the leftmost and rightmost aligned symbols
    UINT32  J;      // bits 0-30: 0-based signed offset from start of R sequence; bit 31: set if R is reverse complement
    INT16   Il;     // leftmost aligned symbol (0-based offset from start of Q sequence)
    INT16   Ir;     // rightmost aligned symbol (0-based offset from start of Q sequence)
    UINT8   subId;  // R-sequence "subunit ID" (e.g. chromosome number)
    UINT8   nTBO;   // number of traceback origins in the same DP band
    INT16   cb;     // number of bytes in the BRLEA binary string
};
#pragma pack(pop)

#pragma pack(push, 2)
struct __align__(8) Qwarp
{
    INT64   sqId[CUDATHREADSPERWARP];   // unique IDs for the sequences in the Qwarp
    INT16   N[CUDATHREADSPERWARP];      // N for each Q sequence in the Qwarp
    UINT32  nAn[CUDATHREADSPERWARP];    // number of nongapped alignments for each Q sequence
    UINT32  nAg[CUDATHREADSPERWARP];    // number of gapped alignments for each Q sequence
    UINT16  nAc[CUDATHREADSPERWARP];    // tentative number of concordant mappings
    INT32   tlen[CUDATHREADSPERWARP];   // TLEN (fragment length) for tentative concordant mappings
    UINT32  maxC[CUDATHREADSPERWARP];   // maximum seed-and-extend seed coverage for each Q sequence
    UINT32  ofsQAI[CUDATHREADSPERWARP]; // offset into QAI buffer
    UINT32  ofsQi;                      // offset (into QiBuffer) of the interleaved Q data for this Qwarp
    INT32   celQmax;                    // maximum number of 64-bit elements in the Q sequences
    UINT32  ofsS;                       // offset (into scoring matrix buffer) of the first scoring matrix for this Qwarp
    INT32   celSmax;                    // maximum number of cells in a scoring matrix for this Qwarp
    INT16   Nmax;                       // maximum N in the Qwarp
    INT16   nQ;                         // number of Q sequences referenced in this Qwarp struct
    INT16   wcgsc;                      // worst-case gap space count for this Qwarp
};
#pragma pack(pop)

/// <summary>
/// Struct <c>CudaGmemLayout</c> defines the layout of CUDA global memory buffers that contain lookup table (hash table) data.
/// </summary>
struct CudaGmemLayout
{
    /* This array of structs describes the "partitions" of a single H or J lookup table that is partitioned
        across multiple GPUs when CUDA "peer access" (i.e., NVlink) is available.

       This entire struct is copied verbatim into CUDA constant memory, so we are obligated to preallocate
        a maximum-size array (one item for each bit in a UINT32).  But who cares about an extra kilobyte or two...
    */
    CudaGmemLayoutType          layoutType;         // memory layout type
    size_t                      cel;                // total number of UINT32 items in the buffer
    size_t                      celPerPartition;    // number of UINT32 items in each partition
    INT32                       nPartitions;        // 1: part[deviceOrdinal]; >1: part[offset/celPerPartition]
    WinGlobalPtrBase*           HostBuf;            // cgltPinned: CudaPinnedPtr; otherwise WinGlobalPtr
    char                        tag[4];             // buffer contents (e.g., "Hn", "Jn", "Hg", "Jg")

    struct CudaGmemPartitionInfo
    {
        CudaGlobalPtr<UINT32>*  PartBuf;    // per-device global memory buffer; NULL if cgltPinned
        UINT32*                 p;          // device pointer (i.e., DeviceBuf->p)
        size_t                  ofs;        // offset into buffer data
        size_t                  cel;        // number of items in buffer partition
    } part[CUDAMAXGPUS];
};
#pragma endregion

class AriocBase
{
    private:
        static const DWORD LUT_LOAD_TIMEOUT = (6*60000);    // 6 minutes
        static const UINT32 celPad = 4*1024;                // number of 64-bit values used as padding for R sequence buffer

    protected:
        INT64                           m_cbTotalAvailable; // total amount of available physical memory reported by the OS
        INT64                           m_cbTotalUsed;      // total bytes of host memory used for all buffers
        INT64                           m_cbTotalPinned;    // total bytes of host memory used for pinned (page-locked) buffers
        WinGlobalPtr<UINT64>            m_Rbuf;             // consolidated host buffer for all R sequence data
        CudaPinnedPtr<UINT64>           m_Rpinned;          // consolidated page-locked buffer for all R sequence data
        ReferenceDataInfo               m_rdiNongapped;     // reference data info for nongapped alignments
        ReferenceDataInfo               m_rdiGapped;        // reference data info for gapped alignments
        HiResTimer                      m_hrt;              // timer
        
    public:
        INT32                           nLPs;           // number of CPU threads (hyperthreads)
        INT16                           nGPUs;          // number of GPU devices
        INT64                           nReadsInFiles;  // expected number of reads (mates) to align
        volatile INT64                  nReadsAligned;  // number of reads (mates) aligned
        UINT32                          RefId;          // R sequence reference ID (AriocE srcId)
        INT32                           minSubId;       // minimum subunit ID (e.g. chromosome number) value
        INT32                           maxSubId;       // maximum subunit ID value
        UINT32                          BatchSize;      // maximum number of Q sequences in one QBatch instance
        WinGlobalPtr<UINT64>            ofsRplus;       // offsets to R+ (forward) sequence data
        WinGlobalPtr<UINT64>            ofsRminus;      // offsets to R- (reverse complement) sequence data
        WinGlobalPtr<UINT32>            M;              // M (sequence length) for each R sequence
        UINT64*                         R;              // consolidated host buffer for all R sequence data
        UINT64                          celR;           // number of 64-bit elements in the R-sequence buffer
        bool                            usePinnedR;     // flag set to place R buffer in CUDA "pinned" memory
        CudaGmemLayout                  cglHn;          // CUDA peer buffer info for nongapped H table
        CudaGmemLayout                  cglJn;          // CUDA peer buffer info for nongapped J table
        CudaGmemLayout                  cglHg;          // CUDA peer buffer info for gapped H table
        CudaGmemLayout                  cglJg;          // CUDA peer buffer info for gapped J table
        InputFileGroup*                 pifgQ;          // Q sequence files
        A21SpacedSeed                   a21ss;          // supports spaced seeds for nongapped alignment
        A21HashedSeed                   a21hs;          // supports hashed seed-and-extend seeds for gapped alignment
        AriocAlignmentScorer            aas;            // alignment score parameters
        
#if TODO_CHOP_IF_UNUSED
        WinGlobalPtr<baseARowWriter*>   SAMwriter;      // for SAM files: one output-file writer instance for each alignment-result reporting flag
        WinGlobalPtr<baseARowWriter*>   SBFwriter;      // for SBF files: one output-file writer instance for each alignment-result reporting flag
        WinGlobalPtr<baseARowWriter*>   TSEwriter;      // for TSE files: one output-file writer instance for each alignment-result reporting flag
        WinGlobalPtr<baseARowWriter*>   KMHwriter;      // for KMH files: one output-file writer instance for each alignment-result reporting flag
#endif

        AriocAppMainBase*               paamb;          // reference to this application's AriocAppMainBase instance
        SAMHDBuilder                    SAMhdb;         // SAM file header builder
        RaiiInverseSemaphore            semLoadLUTs;    // synchronization object for lookup table load
        RaiiInverseSemaphore            semUnloadLUTs;  // synchronization object for lookup table unload
        RaiiInverseSemaphore            semMainLoop;    // synchronization object for joining main-loop CPU threads
        bool                            doMainLoopJoin; // flag to indicate whether to join CPU threads in main loop
        bool                            preferLoadRix;  // true: if possible, use baseLoadRix (optimized interleaved R sequence loader)
        INT32                           TLENbias;       // adjustment to TLEN criterion for concordant mapping

        UINT64                          nTLEN;          // total number of concordant paired-end mappings (i.e., having a valid TLEN)
        UINT64                          sumTLEN;        // sum of valid TLEN values (for computing TLEN standard deviation)
        UINT64                          sosTLEN;        // sum of squares (for computing TLEN standard deviation)
        INT32                           iMeanTLEN;      // mean and standard deviation for TLEN for concordant nongapped mappings
        double                          dMeanTLEN;
        double                          stdevTLEN;
        INT32                           KmerSize;       // kmer size (for kmer hashing)
        CIGARfmtType                    CFT;            // CIGAR format type
        MDfmtType                       MFT;            // MD format type
        UINT32                          StrandsPerSeed; // 2: seeding both Qf and Qrc; 1: seeding only Qf
        static AriocAppMetrics          aam;
        static AA<AriocTaskUnitMetrics> tum;
        tuWatchdog                      Watchdog;

    private:
        void initARowWriters( WinGlobalPtr<OutputFileInfo>& ofi );
        // TODO: CHOP: void zapARowWriters( WinGlobalPtr<baseARowWriter*>& arwList );
        void sniffHostResources( INT32 maxDOP );
        void sniffGpuResources( UINT32 gpuMask );

        INT64 getSqIdFromFile( char* fileSpec );
        INT64 getSqIdAndMfromFile( char* fileSpec );
        INT64 getRfromFile( INT64 ofs, char* fileSpec );
        INT64 appendRpadding( INT64 ofs );
        tuLoadLUT<UINT32>* loadHJhostBuf( CudaGmemLayout& _cgl, char* _fileSpec, bool _serialLUTinit );
        void initHJBufferLayout( CudaGmemLayout& _cgl, tuLoadLUT<UINT32>* _ptuLoadLUT, const char* _tag );
        static int fileSubIdComparer( const void*, const void* );
        void emitCudaGmemLayout( CDPrintFlags _cdpf, CudaGmemLayout* _pcgl );

    protected:
        void loadR( void );
        void loadHJ( void );
        void releaseGpuResources( void );
        void flushARowWriters( WinGlobalPtr<baseARowWriter*>& arwList );

    public:
        AriocBase( const char* _pathR, A21SpacedSeed _a21ss, A21HashedSeed _a21hs, AlignmentControlParameters _acp, AlignmentScoreParameters _asp,
                    InputFileGroup* _pifgQ,
                    WinGlobalPtr<OutputFileInfo>& _ofi,
                    UINT32 _gpuMask, INT32 _maxDOP, UINT32 _batchSize, INT32 _kmerSize, CIGARfmtType _cigarFmtType, MDfmtType _mdFmtType,
                    AriocAppMainBase* paamb );
        virtual ~AriocBase( void );
        void UpdateProgress( UINT32 nQ );
        void DumpCudaGmemLayout( CDPrintFlags _cdpf );

        static AriocTaskUnitMetrics* GetTaskUnitMetrics( const char* _key, const char* _baseKey = NULL );
};
