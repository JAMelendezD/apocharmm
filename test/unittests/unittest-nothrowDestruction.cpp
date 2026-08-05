// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// ENDLICENSE

#include "catch.hpp"

#include "BEDSSubscriber.h"
#include "CompositeSubscriber.h"
#include "CudaBondedForce.h"
#include "CudaContainer.h"
#include "CudaEMap.h"
#include "CudaEnergyVirial.h"
#include "CudaIntegrator.h"
#include "CudaLangevinPistonIntegrator.h"
#include "CudaLangevinThermostatIntegrator.h"
#include "CudaNeighborList.h"
#include "CudaNeighborListBuild.h"
#include "CudaNeighborListSort.h"
#include "CudaP21NeighborListBuild.h"
#include "CudaPMEDirectForce.h"
#include "CudaPMERecip.h"
#include "CudaTopExcl.h"
#include "DcdSubscriber.h"
#include "DeviceVector.h"
#include "DualTopologySubscriber.h"
#include "DynaSubscriber.h"
#include "EDSSubscriber.h"
#include "FEPSubscriber.h"
#include "Force.h"
#include "ForceManager.h"
#include "HarmonicCenterOfMassRestraintForce.h"
#include "HarmonicRestraintForce.h"
#include "Matrix3d.h"
#include "RestartSubscriber.h"
#include "StateSubscriber.h"
#include "Subscriber.h"
#include "XYZSubscriber.h"

#include <type_traits>
#include <utility>

static_assert(std::is_nothrow_destructible<DeviceVector<float>>::value);
static_assert(std::is_nothrow_destructible<CudaContainer<float>>::value);
static_assert(std::is_nothrow_destructible<Force<long long int>>::value);
static_assert(std::is_nothrow_destructible<Matrix3d<float>>::value);
static_assert(std::is_nothrow_destructible<CudaEnergyVirial>::value);

static_assert(std::is_nothrow_destructible<
              HarmonicRestraintForce<long long int, float>>::value);
static_assert(std::is_nothrow_destructible<
              HarmonicCenterOfMassRestraintForce<long long int, float>>::value);

static_assert(noexcept(std::declval<CudaIntegrator &>().~CudaIntegrator()));
static_assert(
    std::is_nothrow_destructible<CudaLangevinPistonIntegrator>::value);
static_assert(
    std::is_nothrow_destructible<CudaLangevinThermostatIntegrator>::value);

static_assert(std::is_nothrow_destructible<ForceManager>::value);
static_assert(std::is_nothrow_destructible<CudaEMap>::value);
static_assert(std::is_nothrow_destructible<CudaTopExcl>::value);
static_assert(
    std::is_nothrow_destructible<CudaBondedForce<long long int, float>>::value);
static_assert(std::is_nothrow_destructible<
              CudaPMEDirectForce<long long int, float>>::value);

static_assert(std::is_nothrow_destructible<CudaNeighborListBuild<32>>::value);
static_assert(std::is_nothrow_destructible<CudaP21NeighborListBuild>::value);
static_assert(std::is_nothrow_destructible<CudaNeighborListSort>::value);
static_assert(std::is_nothrow_destructible<CudaNeighborList<32>>::value);

static_assert(
    std::is_nothrow_destructible<CudaPMERecip<int, float, float2>>::value);
static_assert(std::is_nothrow_move_constructible<
              CudaPMERecip<int, float, float2>>::value);

static_assert(noexcept(std::declval<Subscriber &>().~Subscriber()));
static_assert(std::is_nothrow_destructible<DcdSubscriber>::value);
static_assert(std::is_nothrow_destructible<DynaSubscriber>::value);
static_assert(std::is_nothrow_destructible<BEDSSubscriber>::value);
static_assert(std::is_nothrow_destructible<EDSSubscriber>::value);
static_assert(std::is_nothrow_destructible<FEPSubscriber>::value);
static_assert(std::is_nothrow_destructible<RestartSubscriber>::value);
static_assert(std::is_nothrow_destructible<StateSubscriber>::value);
static_assert(std::is_nothrow_destructible<XYZSubscriber>::value);
static_assert(std::is_nothrow_destructible<DualTopologySubscriber>::value);
static_assert(std::is_nothrow_destructible<CompositeSubscriber>::value);

TEST_CASE("CUDA resource owners have nothrow destructors") { SUCCEED(); }
