/*
  AriocMetricsBase.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#pragma once
#define __AriocMetricsBase__

/// <summary>
/// Class <c>AriocMetricsBase</c> implements base functionality for derived performance-metrics class implementations.
/// </summary>
class AriocMetricsBase
{
    protected:
        static const double m_msPerSec;
        static const double m_usPerSec;
        static const double m_bytesPerMB;
        static const double m_bytesPerGB;

    protected:
        AriocMetricsBase( void );
        virtual ~AriocMetricsBase( void );
        UINT64 usToMs( const UINT64 us );
};
