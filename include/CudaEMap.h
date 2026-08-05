// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Andrew Simmonett, Samarjeet Prasad
//
// ENDLICENSE

#pragma once
#include "CharmmContext.h"
#include "CudaContainer.h"
#include <cuda_runtime.h>
#include <memory>

class CudaEMap {
public:
  // for now we are assuming that the grid is cubic
  // and set to 100x100x100.
  // this can be changed later
  CudaEMap(std::shared_ptr<CharmmContext> context)
      : context(context), stream(nullptr), nx(100), ny(100), nz(100) {
    cudaCheck(cudaStreamCreate(&stream));

    try {
      atomicMasses = context->getForceManager()->getPsf()->getMasses();
    } catch (...) {
      destroy_cuda_stream_noexcept(&stream);
      throw;
    }
  }

  ~CudaEMap(void) noexcept { destroy_cuda_stream_noexcept(&stream); }

  void generate();
  void dealloc(void);

private:
  std::shared_ptr<CharmmContext> context;

  CudaContainer<double> atomicMasses;

  // cudastream
  cudaStream_t stream;
  // fft
  // fftx, ffty, fftz

  int nx, ny, nz;
};
