/*
  CpuInfo.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma region constructors/destructor
// default constructor
CpuInfo::CpuInfo() : Description()
{
    /* Look for the CPU info string maintained by the OS.  If an error occurs, do nothing. */

#ifdef _WIN32
    HKEY hKey;

    // open the registry key
    long rval = RegOpenKeyEx( HKEY_LOCAL_MACHINE,
        "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0",
        0,
        KEY_READ,
        &hKey );
    if( rval != ERROR_SUCCESS )
        return;

    // query the registry value
    DWORD cb;
    rval = RegQueryValueEx( hKey, "ProcessorNameString", NULL, NULL, reinterpret_cast<LPBYTE>(this->Description), &cb );
    if( rval != ERROR_SUCCESS )
        return;
#endif

#ifdef __GNUC__
    RaiiFile ProcCpuInfo;

    const char* filespec = "/proc/cpuinfo";
    ProcCpuInfo.OpenReadOnly( const_cast<char*>(filespec) );
    if( ProcCpuInfo.Handle < 0 )
        return;

    /* /proc/cpuinfo starts like this:
    
        processor       : 0
        vendor_id       : GenuineIntel
        cpu family      : 6
        model           : 63
        model name      : Intel(R) Xeon(R) CPU E5-2695 v3 @ 2.30GHz

       We just stupidly look for "model name" and then grab whatever lies between the colon and the newline.
    */
    char buf[256];
    INT64 cbRead = ProcCpuInfo.Read( buf, sizeof buf );
    if( cbRead != (sizeof buf) )
        return;

    char* p = strstr( buf, "model name" );
    if( p == NULL )
        return;

    p = strchr( p, ':' );
    if( p == NULL )
        return;

    char* pEOL = strpbrk( ++p, "\r\n" );
    if( pEOL == NULL )
        return;

    *pEOL = '\0';

    strcpy_s( this->Description, sizeof this->Description, p );
#endif

    UnduplicateSpaces( this->Description );
}

// destructor
CpuInfo::~CpuInfo()
{
}
#pragma endregion

