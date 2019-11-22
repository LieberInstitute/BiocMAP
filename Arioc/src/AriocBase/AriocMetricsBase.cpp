/*
   AriocMetricsBase.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

#pragma region static member variable definitions
const double AriocMetricsBase::m_msPerSec = 1000;
const double AriocMetricsBase::m_usPerSec = 1000000;
const double AriocMetricsBase::m_bytesPerMB = 1024*1024;
const double AriocMetricsBase::m_bytesPerGB = 1024*1024*1024;
#pragma endregion

#pragma region constructor/destructor
/// default constructor
AriocMetricsBase::AriocMetricsBase()
{
}

/// destructor
AriocMetricsBase::~AriocMetricsBase()
{
}
#pragma endregion

#pragma region protected methods
/// [protected] method usToMs
UINT64 AriocMetricsBase::usToMs( const UINT64 us )
{
    // convert microseconds to milliseconds with rounding to the nearest millisecond
    return (us + 500) / 1000;
}
#pragma endregion
