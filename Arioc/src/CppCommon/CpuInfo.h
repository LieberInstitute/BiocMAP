/*
  CpuInfo.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/

#pragma once
#define __CpuInfo__

/// <summary>
/// Class <c>CpuInfo</c> obtains CPU name and description from the OS.
/// </summary>
class CpuInfo
{
    public:
        char Description[128];

    public:
        CpuInfo( void );
        virtual ~CpuInfo( void );
};
