/*
  QBatchPool.h

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#pragma once
#define __QBatchPool__

/// <summary>
/// Class <c>QBatchPool</c> maintains a pool of QBatch instances.
/// </summary>
class QBatchPool
{
    static const UINT32 POOL_TIMEOUT = 5 * 60000;   // 5 minutes

    private:
        WinGlobalPtr<QBatch*>   m_qb;
        RaiiSemaphore           m_semPool;
        RaiiSyncEventObject     m_evPoolFull;

    public:
        INT16                   Nmax;   // maximum read length

    private:
        QBatchPool( void );

    public:
        QBatchPool( INT32 celPool, AriocBase* _pab, GpuInfo* _pgi, INT16 _Nmax );
        ~QBatchPool( void );
        void Enpool( QBatch* pqb );
        QBatch* Depool( InputFileGroup::FileInfo* _pfi );
        bool Wait( UINT32 msTimeout );
};
