/*
  tuTailU.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#pragma once
#define __tuTailU__

#ifndef __tuClassifyU__
#include "tuClassifyU.h"
#endif

/// <summary>
/// Class <c>tuTailU</c> classifies and writes unpaired alignment results to disk files
/// </summary>
class tuTailU : public tuTail
{
    private:
        tuTailU( void );

    protected:
        virtual void classifyAlignmentResults( void );

    public:
        tuTailU( QBatch* pqb );
        virtual ~tuTailU( void );

        virtual baseWriteA* instantiateARowFormatter( OutputFormatType _oft, INT32 _iBuf, bool _wantCounts );
};

