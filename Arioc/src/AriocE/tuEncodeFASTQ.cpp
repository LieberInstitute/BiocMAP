/*
  tuEncodeFASTQ.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma region constructors and destructor
/// <summary>
/// Encodes FASTQ-formatted sequence data.
/// </summary>
/// <param name="psip">Reference to a common parameter structure</param>
/// <param name="sqCat">sequence category (+ or - strand)</param>
/// <param name="iInputFile">index of input file</param>
/// <param name="psem">Reference to an <c>RaiiSemaphore</c> instance</param>
/// <param name="pscw">Reference to a <c>SAMConfigWriter</c> instance (may be NULL)</param>
tuEncodeFASTQ::tuEncodeFASTQ( AriocEncoderParams* psip, SqCategory sqCat, INT16 iInputFile, RaiiSemaphore* psem, SAMConfigWriter* pscw ) : baseEncode(psip,sqCat,iInputFile,psem,pscw)
{
}

/// [public] destructor
tuEncodeFASTQ::~tuEncodeFASTQ()
{
}
#pragma endregion

#pragma region private methods
/// [private] method summarizeReadLengthDistribution
void tuEncodeFASTQ::summarizeReadLengthDistribution()
{
    // do nothing if no reads were processed
    if( m_nN.n == 0 )
        return;

    // minimum read length
    UINT32 minN = 0;
    UINT32 N = 0;
    for( ; N<=m_nN.n; ++N )
    {
        if( m_nN.p[N] )
        {
            minN = N;
            break;
        }
    }

    // maximum read length
    UINT32 maxN = 0;
    N = m_nN.n;
    do
    {
        if( m_nN.p[N] )
        {
            maxN = N;
            break;
        }
    }
    while( --N != 0 );

    // read count, maximum count, mean read length
    UINT32 nN = 0;
    UINT32 maxNcount = 0;
    UINT64 sumN = 0;
    for( UINT32 N=minN; N<=maxN; ++N )
    {
        nN += m_nN.p[N];
        maxNcount = max2(maxNcount, m_nN.p[N]);
        sumN += m_nN.p[N] * static_cast<UINT64>(N);
    }

    double avgN = static_cast<double>(sumN) / nN;

    // standard deviation
    double ss = 0;
    for( UINT32 N=minN; N<=maxN; ++N )
        ss += m_nN.p[N] * ((N-avgN) * (N-avgN));
    double sdN = sqrt( ss/nN );

    // emit results
    {
        RaiiCriticalSection<tuEncodeFASTQ> rcs;

        INT32 cchN = static_cast<INT32>(log10( static_cast<double>(maxN) )) + 1;            // number of digits in maximum read length
        INT32 cchNcount = static_cast<INT32>(log10( static_cast<double>(maxNcount) )) + 1;  // number of digits in maximum count

        CDPrint( cdpCDg, "%s: Read-length distribution for %s:", __FUNCTION__, m_inFile.FileSpec.p );
        CDPrint( cdpCDg, "%s: %-*s %s", __FUNCTION__, cchN, "N", "nReads" );
        for( UINT32 N=minN; N<=maxN; ++N )
            CDPrint( cdpCDg, "%s: %*d %*d", __FUNCTION__, cchN, N, cchNcount, m_nN.p[N] );

        CDPrint( cdpCD1, "%s: %s: %u reads [%u-%u mean (sd) %2.1f (%2.1f)]", __FUNCTION__, m_inFile.FileSpec.p, nN, minN, maxN, avgN, sdN );
    }

    // save the maximum read length
    m_psip->maxN = maxN;
}

/// [private] method computeBaseTripletCount
void tuEncodeFASTQ::computeBaseTripletCount( char* _p, INT32 _cch )
{
    const char* const alphabet = "ACGT";
    UINT64 bits = { 0 };

    // start with the last three characters in the string
    char* p = _p + (_cch - 1);
    const char* ip = strchr( alphabet, toupper( *p ) );
    INT32 i2 = ip ? static_cast<INT32>(ip - alphabet) : -1;

    --p;
    ip = strchr( alphabet, toupper( *p ) );
    INT32 i1 = ip ? static_cast<INT32>(ip - alphabet) : -1;

    while( --p >= _p )
    {
        ip = strchr( alphabet, toupper( *p ) );
        INT32 i0 = ip ? static_cast<INT32>(ip - alphabet) : -1;

        // compute the bit index; this value will be negative if any of the 3 characters was not in the alphabet
        INT32 shl = (i2 << 4) | (i1 << 2) | i0;

        // set the bit corresponding to the triplet sequence
        if( shl >= 0 )
            bits |= static_cast<UINT64>(1) << shl;

        // iterate
        i2 = i1;
        i1 = i0;
    };

#if TODO_MAYBE_SOMEDAY
    for( INT32 i=(_cch-3); i>=0; --i )
    {
        // convert the three-letter substring at the i'th position in p into a 6-bit value
        INT32 ix = strchr( "ACGT", *


            /*
        A:  41h     0100 0001
        C:  43h     0100 0011
        G:  47h     0100 0111
        T:  54h     0101 0100


            N:  4Eh 0100 1110

            */

            /* Get 3 bytes (3 ASCII characters) as a 32-bit integer.  This type punning actually delivers 4 bytes,
                so theoretically the read at the end of the character string "overflows" the buffer -- but we know
                that the buffer is always bigger than the sequence characters because there's more FASTQ data
                immediately afterward.

               This gives us a 3-character triplet in the 24 low-order bits of the integer.
            */
            INT32 c3 = *reinterpret_cast<INT32*>(_p+i);

        /* Isolate the three characters by zeroing irrelevant bits:
            bits 24..31:  character following the three characters of interest
            bit  21:      lowercase, 3rd character
            bit  13:      lowercase, 2nd character
            bit   5:      lowercase, 1st character
        */
        c3 &= 0x00DFDFDF;

        /* Subtract 'A' from each character.

           At this point the high-order byte must be zero and each of the 3 low-order bytes should contain
            one of the following values:
                                          shr 1      add bit 3 shr 3
                 0 (A)    00000000        00000000   00000000
                 3 (C)    00000011        00000001   00000001
                 6 (G)    00000110        00000011   00000011
                19 (T)    00010011        00001001   00001010
                13 (N)    00011110        00001111   00010000

                ????????
           No error checking here!
        */

        }
#endif

    // save the number of different triplets
    m_btcCurrent = static_cast<UINT8>(__popcnt64( bits ));
}
#pragma endregion

#pragma region virtual method implementations
/// <summary>
/// Encodes FASTQ-formatted sequence data.
/// </summary>
void tuEncodeFASTQ::main()
{
    CDPrint( cdpCD3, "%s...", __FUNCTION__ );

    /* Set up to write a .CFG file for the input sequence file.  If the caller does not provide a
        SAMConfigWriter instance, we create one here.
    */
    if( m_pscw == NULL )
        m_pscw = new SAMConfigWriter( m_psip, m_psip->OutFilespecStubSq, m_baseName );

    // assign a 0-based ordinal to each sequence
    INT64 readId = 0;
    UINT64 nSqIn = 0;
    UINT64 nSqEncoded = 0;

    // read the input file in chunks
    m_inFile.Seek( 0, SEEK_SET );
    INT64 cbRead = readInputFile( m_inBuf.p, m_inBuf.cb );

    char* pLimit = m_inBuf.p + cbRead;
    char* pEOL = m_inBuf.p;
    char* p = findSOL( pEOL, &pLimit );

    while( p && (p < pLimit) )
    {
        // find the end of the next line of data and append a null terminator
        char* pEOL = findEOL( &p, &pLimit );

        /* compute the hash threshold for encoding the current sequence; that is, we only encode the current sequence if
            the 32-bit hash of its readId is below a threshold determined by the user-specified sampling ratio */
        bool encodeSq = (AriocE::Hash6432( readId ) < m_psip->EncodeThreshold);

        // copy the sequence metadata into its output buffer
        p++ ;                                           // point past the '@' symbol

        if( encodeSq )
            (this->*m_writeSqm)(readId, p, static_cast<INT32>(pEOL-p));       // copy the data

        // copy and encode the sequence data
        p = findSOL( pEOL+1, &pLimit );
        INT32 cchSq = 0;                                // count the sequence symbols
        while( p && (*p != '+') )
        {
            // find the end of the next line of data and append a null terminator
            pEOL = findEOL( &p, &pLimit );

            // count sequence symbols
            INT32 cch = static_cast<INT32>(pEOL-p);
            cchSq += cch;

            // copy the raw and encoded sequence data into their output buffers
			if( encodeSq )
			{
				(this->*m_writeRaw)(readId, p, cch);
				computeBaseTripletCount( p, cch );
			}

            // advance to the next line in the input file
            p = findSOL( pEOL+1, &pLimit );
        }

        // sanity check
        if( p == NULL )
            throw new ApplicationException( __FILE__, __LINE__, "unexpected end of file: %s", m_inFile.FileSpec.p );

        // skip past the '+' row
        pEOL = findEOL( &p, &pLimit );
        p = findSOL( pEOL+1, &pLimit );
        INT32 cchQ = 0;                                // count the quality scores

        // process the quality scores
        bool readingQualityScores = true;
        do
        {
            // we can't simply look for '@' here because this symbol may appear as the first character in a line of quality scores
            while( p && (cchQ < cchSq) )
            {
                // find the end of the next line of data and append a null terminator
                pEOL = findEOL( &p, &pLimit );

                // count quality scores
                INT32 cch = static_cast<INT32>(pEOL-p);
                cchQ += cch;

                // copy the quality scores into their output buffer
                if( encodeSq )
                    (this->*m_writeSqq)(readId, p, cch);

                // advance to the next line in the input file
                p = findSOL( pEOL+1, &pLimit );
            }

            if( (p == NULL) || (cchQ >= cchSq) )
                readingQualityScores = false;
            else
            {
                pEOL = findEOL( &p, &pLimit );
                cchQ += static_cast<INT32>(pEOL-p);                // count quality scores
                p = findSOL( pEOL+1, &pLimit );
            }
        }            
        while( readingQualityScores );

        // terminate the row data in the raw and quality-score output buffers
        if( encodeSq )
        {
            (this->*m_endRowRaw)();
            endRowSqqQ();
        }

        // increment the counts
        readId++;           // ordinal (0-based)
        nSqIn++;            // number of sequences read
        if( encodeSq )
            nSqEncoded++ ;  // number of sequences encoded

        // accumulate the read-length distribution
        if( static_cast<UINT32>(cchSq) > m_nN.n )
        {
            m_nN.Realloc( cchSq+1, true );
            m_nN.n = static_cast<UINT32>(cchSq);
        }
        m_nN.p[cchSq]++ ;
    }

    // flush the output buffers
    flushSqm();
    flushRaw( true );
    flushSqq( true );
    flushA21();
    flushKmers();

    // performance metrics
    InterlockedExchangeAdd( &m_psip->nSqIn, nSqIn );
    InterlockedExchangeAdd( &m_psip->nSqEncoded, nSqEncoded );

    // summarize the read-length distribution
    summarizeReadLengthDistribution();

    // write the .cfg file (query sequences only)
    if( m_psip->Nenc == NencodingQ )
        writeConfigFile();

    // signal that this thread has completed its work
    m_psemComplete->Release( 1 );

    CDPrint( cdpCD3, "%s completed", __FUNCTION__ );
}
#pragma endregion
