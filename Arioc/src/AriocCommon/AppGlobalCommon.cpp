/*
  AppGlobalCommon.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#ifndef __AriocVersion_Includes__
#include "../AriocVersion/AriocVersion.h"
#endif

#pragma region static variable initialization
char AppGlobalCommon::MachineName[MAX_COMPUTERNAME_LENGTH+2];
bool AppGlobalCommon::WantEnterAtExit = false;
char AppGlobalCommon::m_dtBuild[20] = { 0 };
char AppGlobalCommon::m_fmtDefaultCfgFile[] =
#ifdef _WIN32
                                  "%s.cfg";
#endif
#ifdef __GNUC__
                                  "./%s.cfg";
#endif
#pragma endregion

#pragma region constructor/destructor
/// [protected] default constructor
AppGlobalCommon::AppGlobalCommon() : m_pamc(NULL)
{
    // do not use this constructor
    throw new ApplicationException( __FILE__, __LINE__, "(not implemented)" );
}

/// constructor (AppMainCommon*)
AppGlobalCommon::AppGlobalCommon( AppMainCommon* pamc ) : m_pamc(pamc),
                                                          m_szExeFilespec(), m_szAppName(), m_szAppVersion(), m_szAppCfg(), m_szAppLegal()
{
}

/// destructor
AppGlobalCommon::~AppGlobalCommon()
{
}
#pragma endregion

#pragma region protected methods
/// [protected] method getStartupBannerDetails
void AppGlobalCommon::getStartupBannerDetails( const char* _argv0 )
{
    /* The executable file specification, application name, and application version are handled a bit differently
        in Windows and in Linux:
        - in Windows, the buffers get updated here with strings extracted from the "fixed file info"
            in the executable file
        - in Linux, the app version is copied from an include file that gets updated in the Windows Visual Studio
            project (i.e., for distribution builds, the Windows version should be compiled before the Linux
            version)
    */

#if defined(_WIN32)
    // extract the file version info
    if( 0 == GetModuleFileName( NULL, m_szExeFilespec, FILENAME_MAX ) )
    {
        DWORD dwErr = GetLastError();
        throw new ApplicationException( __FILE__, __LINE__, "GetModuleFileName returned 0 (GetLastError returned %u)", dwErr );
    }

    DWORD dwHandle;
    UINT32 cb = GetFileVersionInfoSize( m_szExeFilespec, &dwHandle );
    WinGlobalPtr<INT8> vi( cb, true );
    GetFileVersionInfo( m_szExeFilespec, dwHandle, cb, vi.p );
    if( cb == 0 )
        throw new ApplicationException( __FILE__, __LINE__, "missing file version info" );

    // the application name is the "internal name" in the Windows version info
    char* p;
    VerQueryValue( vi.p, "\\StringFileInfo\\040904E4\\InternalName", reinterpret_cast<void**>(&p), &cb );
    strcpy_s( m_szAppName, APPNAME_BUFSIZE, p );

    // the major and minor version numbers are in the "product version" string
    VerQueryValue( vi.p, "\\StringFileInfo\\040904E4\\FileVersion", reinterpret_cast<void**>(&p), &cb );
    strcpy_s( m_szAppVersion, APPVER_BUFSIZE, p );

    // the copyright info is in the "legal copyright" string
    VerQueryValue( vi.p, "\\StringFileInfo\\040904E4\\LegalCopyright", reinterpret_cast<void**>(&p), &cb );
    strcpy_s( m_szAppLegal, APPLEGAL_BUFSIZE, p );
#endif

#ifdef __GNUC__
    strcpy_s( m_szAppVersion, sizeof m_szAppVersion, VER_FILEVERSION_STR );
    strcpy_s( m_szAppLegal, sizeof m_szAppLegal, VER_LEGALCOPYRIGHT_STR );
    realpath( _argv0, m_szExeFilespec );
#endif
}

/// [protected virtual] method parseCommandTail
void AppGlobalCommon::parseCommandTail( int argc, char* argv[] )
{
    /* Parse the command line:    
        argv[0]: executable filename
        argv[1]: configuration filename
    */

    switch( argc )
    {
        case 1:
            sprintf_s( m_szAppCfg, sizeof m_szAppCfg, m_fmtDefaultCfgFile, m_szAppName );
            break;

        case 2:
            strcpy_s( m_szAppCfg, sizeof m_szAppCfg, argv[1] );
            break;

        default:
            throw new ApplicationException( __FILE__, __LINE__, "invalid command-line arguments (syntax: %s [configuration_filename])", m_szAppName );
    }
}
#pragma endregion

#pragma region public methods
/// [public] method Run
int AppGlobalCommon::Run( int argc, char* argv[], const char* defaultAppName )
{
    int rval = 0;

    UINT32 cbMachineName = sizeof AppGlobalCommon::MachineName;
    GetComputerName( AppGlobalCommon::MachineName, reinterpret_cast<LPDWORD>(&cbMachineName) );

    getStartupBannerDetails( argv[0] );
    if( m_szAppName[0] == 0 )
        strcpy_s( m_szAppName, sizeof m_szAppName, defaultAppName );

    CDPrint( cdpCD0, "%s v%s (%s) +%s"
                   , m_szAppName
                   , m_szAppVersion
#ifdef _DEBUG
                   , "debug"
#else
                   , "release"
#endif
                   , VER_REVISIONDATETIME_STR      
           );

    try
    {
        CDPrint( cdpCD0, m_szAppLegal );
        CDPrint( cdpCD0, " computer name        : %s", AppGlobalCommon::MachineName );
        CDPrint( cdpCD0, " executable file      : %s", m_szExeFilespec );
        if( *m_dtBuild )
            CDPrint( cdpCD0, " compiled             : %s", m_dtBuild );
        CDPrint( cdpCD0, " data type sizes      : int=%llu long=%llu *=%llu Jvalue5=%llu Jvalue8=%llu JtableHeader=%llu", sizeof(int), sizeof(long), sizeof(void*), sizeof(Jvalue5), sizeof(Jvalue8), sizeof(JtableHeader) );

        // parse command-line arguments
        parseCommandTail( argc, argv );

        // verify the existence of the configuration file
        char cfgFileSpec[FILENAME_MAX];
        _fullpath( cfgFileSpec, m_szAppCfg, sizeof cfgFileSpec );

        if( !RaiiFile::Exists( cfgFileSpec ) )
            throw new ApplicationException( __FILE__, __LINE__, "cannot open configuration file '%s'", cfgFileSpec );

        // launch the application
        m_pamc->Init( m_szAppName, m_szAppVersion, cfgFileSpec );
        m_pamc->LoadConfig();   // load the XML-formatted config file
        m_pamc->Launch();       // parse the config file and launch the program
    }
    catch( ApplicationException* pex )
    {
        rval = pex->Dump();
    }

    // conditionally wait for the user to press the Enter key
    CDPrintFilter = cdpCD0;
    if( AppGlobalCommon::WantEnterAtExit )
    {
        CDPrint( static_cast<CDPrintFlags>(cdpConsole|0x00000001), "Done.  Press Enter to halt..." );
        getchar();
    }
   
    CDPrint( cdpCD0, "%s ends (%d)", m_szAppName, rval );
    return rval;
}

/// [public static] method SaveDtBuild
bool AppGlobalCommon::SaveBuildTimestamp( const char* _dtBuild )
{
    /* The caller should be the initializer for a C-style static variable in the compilation unit
        whose timestamp is interesting (e.g. AriocU.cpp or AriocP.cpp).  Example:

            static bool isTimestampInitialized = AppGlobalCommon::SaveBuildTimestamp( __TIMESTAMP__ );

       This implementation then gets called during static initialization time, i.e., before the program
        starts running and AppGlobalCommon::Run() gets called.
    */

    // make a local copy of the specified string
    char buf[64];
    strcpy_s( buf, sizeof buf, _dtBuild );

    /* Split the string.
    
       The string comes from the C++ compiler macro __TIMESTAMP__ which is formatted as Ddd Mmm Date hh::mm::ss yyyy, e.g.

            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
            T  u  e     N  o  v        5    1  2  :  1  8  :  3  5     2  0  1  9  \0
    */
    buf[7] = buf[10] = buf[19] = 0;

    // convert the month name to its corresponding number
    const char* MM = (0 == strcmp( buf+4, "Jan" )) ? "01" :
                     (0 == strcmp( buf+4, "Feb" )) ? "02" :
                     (0 == strcmp( buf+4, "Mar" )) ? "03" :
                     (0 == strcmp( buf+4, "Apr" )) ? "04" :
                     (0 == strcmp( buf+4, "May" )) ? "05" :
                     (0 == strcmp( buf+4, "Jun" )) ? "06" :
                     (0 == strcmp( buf+4, "Jul" )) ? "07" :
                     (0 == strcmp( buf+4, "Aug" )) ? "08" :
                     (0 == strcmp( buf+4, "Sep" )) ? "09" :
                     (0 == strcmp( buf+4, "Oct" )) ? "10" :
                     (0 == strcmp( buf+4, "Nov" )) ? "11" :
                     (0 == strcmp( buf+4, "Dec" )) ? "12" :
                     "\0\0";

    // if the month name is unrecognized, use the timestamp string verbatim
    if( *MM == 0 )
    {
        strcpy_s( m_dtBuild, sizeof m_dtBuild, _dtBuild );
        return false;
    }

    /* Build an ISO8601-formatted (yyyy-MM-ddTHH:mm:ss) copy of the timestamp

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
        T  u  e     N  o  v        5    1  2  :  1  8  :  3  5     2  0  1  9  \0

       becomes
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
        2  0  1  9  -  1  1  -  0  5 T  1  2  :  1  8  :  3  5  \0

    */
    
    // replace leading spaces with zeroes
    if( buf[ 8] == ' ' ) buf[ 8] = '0';
    if( buf[11] == ' ' ) buf[11] = '0';
    if( buf[14] == ' ' ) buf[14] = '0';
    if( buf[17] == ' ' ) buf[17] = '0';

    // copy the components of the string
    memcpy_s( m_dtBuild, 20, buf, 20 );
    m_dtBuild[4] = m_dtBuild[7] = '-';
    m_dtBuild[10] = 'T';
    memcpy_s( m_dtBuild+0, 4, buf+20, 4 );
    memcpy_s( m_dtBuild+5, 2, MM, 2 );

    return true;
}
#pragma endregion
