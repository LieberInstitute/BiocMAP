/*
  OutputFileBinding.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.

  Notes:
   The idea here is that each calling thread (i.e., one per GPU) fills its own output buffer; the ARowWriter implementation
    flushes that buffer to disk in a thread-safe manner.
*/
#pragma once
#define __OutputFileBinding__

// forward declaration
class AriocBase;

/// <summary>
/// Class <c>OutputFileBinding</c> binds buffers to output files.
/// </summary>
class OutputFileBinding
{
    private:
        static const INT64 OUTPUT_BUFFER_SIZE = 10 * 1024*1024;     // 10MB

        static WinGlobalPtr<OutputFileBinding*> m_ofb;              // all OutputFileBinding instances
        static AriocBase*                       m_pab;              // reference to this app's one-and-only AriocBase instance
        static INT32                            m_nBuffersPerGPU;

        AlignmentResultFlags                m_arf;
        WinGlobalPtr<WinGlobalPtr<char>*>   m_buf;                  // buffers associated with the output file (one buffer per thread)

    public:
        static INT32 BuffersPerBinding;
        static INT32 ThreadsPerGPU;

        baseARowWriter* pARW;

    private:
        OutputFileBinding();
        OutputFileBinding( AriocBase* _pab, AlignmentResultFlags _arf, baseARowWriter* _parw );
        static void computeOutputBufferCounts( void );

    public:
        virtual ~OutputFileBinding( void );
        // TODO: CHOP IF UNUSED:  void CopyToOutputBuffer( INT32 _iBuf, WinGlobalPtr<char>& _rowBuf );
        char* ReserveOutputBuffer( INT32 _iBuf, UINT32 _cb );
        void CommitOutputBuffer( INT32 _iBuf, UINT32 _cb );

        static void CreateBinding( AriocBase* _pab, AlignmentResultFlags _arf, baseARowWriter* _parw );
        static OutputFileBinding* GetBinding( OutputFormatType _oft, AlignmentResultFlags _arf );
        static void InitOutputBuffers( void );
        static INT32 GetOutputFormatTypes( void );
        static OutputFormatType GetOutputFormatTypeForCounts( void );
        static AlignmentResultFlags GetARFforOFT( OutputFormatType _oft );
        static INT32 GetBufferIndexBase( INT16 _gpuDeviceOrdinal );
        static void Reset( void );
};
