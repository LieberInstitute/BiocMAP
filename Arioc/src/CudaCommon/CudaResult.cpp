/*
  CudaResult.cpp

    Copyright (c) 2015-2019 Johns Hopkins University.  All rights reserved.

    This file is part of the Arioc software distribution.  It is subject to the license terms
    in the LICENSE.txt file found in the top-level directory of the Arioc software distribution.
    The contents of this file, in whole or in part, may only be copied, modified, propagated, or
    redistributed in accordance with the license terms contained in LICENSE.txt.
*/
#include "stdafx.h"

// default constructor
CudaResult::CudaResult() : FileName(), LineNumber(0)
{
}

// constructor (char*)
CudaResult::CudaResult( const char* fileName ) : LineNumber(0)
{
    // save the parameter value
    strcpy_s( this->FileName, sizeof this->FileName, fileName );
}

// destructor
CudaResult::~CudaResult()
{
}

// overload operator=
void CudaResult::operator=( const cudaError_t rval )
{
    // if no error, return immediately
    if( rval == cudaSuccess )
        return;

    // error; throw an exception
    throw new ApplicationException( this->FileName, this->LineNumber, "CUDA runtime API error %d: %s\r\n", rval, cudaGetErrorString(rval) );
}

// method SetDebugInfo
CudaResult& CudaResult::SetDebugInfo( const int lineNumber )
{
    // save the parameter value
    this->LineNumber = lineNumber;

    // return a refererence to this CudaResult instance
    return *this;
}
