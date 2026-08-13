// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Andrew Simmonett, Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include "DeviceVector.h"

#include <cstddef>
#include <vector>
#include <vector_types.h>

template <typename T> class CudaContainer {
public: // Member functions
  CudaContainer(void);
  CudaContainer(const std::size_t count);
  CudaContainer(const std::vector<T> &other);
  CudaContainer(const std::vector<T> &&other);
  CudaContainer(const DeviceVector<T> &other);
  CudaContainer(const DeviceVector<T> &&other);
  CudaContainer(const CudaContainer<T> &other);
  CudaContainer(const CudaContainer<T> &&other);
  ~CudaContainer(void) noexcept = default;

  CudaContainer<T> &operator=(const std::vector<T> &other);
  CudaContainer<T> &operator=(const std::vector<T> &&other);
  CudaContainer<T> &operator=(const DeviceVector<T> &other);
  CudaContainer<T> &operator=(const DeviceVector<T> &&other);
  CudaContainer<T> &operator=(const CudaContainer<T> &other);
  CudaContainer<T> &operator=(const CudaContainer<T> &&other);

public: // Element access
  const T &at(const std::size_t pos) const;
  T &at(const std::size_t pos);

  const T &operator[](const std::size_t pos) const;
  T &operator[](const std::size_t pos);

  const std::vector<T> &getHostArray(void) const;
  std::vector<T> &getHostArray(void);

  const DeviceVector<T> &getDeviceArray(void) const;
  DeviceVector<T> &getDeviceArray(void);

public: // Capacity
  std::size_t size(void) const;
  void shrink_to_fit(void);

public: // Modifiers
  void clear(void);
  void push_back(const T &value);
  void resize(const std::size_t count);
  void set(const std::vector<T> &values);
  void set(const DeviceVector<T> &values);
  void set(const T value);
  void setToValue(const T value);

  void transferToDevice(void);
  void transferToHost(void);
  void transferFromDevice(void);
  void transferFromHost(void);

public:
  void printDeviceArray(void) const;

  // TODO : enable a range-based iterator
private:
  std::vector<T> m_HostArray;
  DeviceVector<T> m_DeviceArray;
};

template class CudaContainer<int>;
template class CudaContainer<int2>;
template class CudaContainer<int3>;
template class CudaContainer<int4>;
template class CudaContainer<unsigned int>;
template class CudaContainer<float>;
template class CudaContainer<float2>;
template class CudaContainer<float3>;
template class CudaContainer<float4>;
template class CudaContainer<long long int>;
template class CudaContainer<longlong2>;
template class CudaContainer<longlong3>;
template class CudaContainer<longlong4>;
template class CudaContainer<unsigned long long int>;
template class CudaContainer<std::size_t>;
template class CudaContainer<double>;
template class CudaContainer<double2>;
template class CudaContainer<double3>;
template class CudaContainer<double4>;
