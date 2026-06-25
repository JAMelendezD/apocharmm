// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "DeviceVector.h"
#include "apo_test_helpers.h"
#include "catch.hpp"

TEST_CASE("ConstructionDestruction") {
  SECTION("DefaultConstructor") {
    DeviceVector<int> v;
    CHECK(v.empty() == true);
    CHECK(v.size() == 0);
    CHECK(v.capacity() == 0);
    CHECK(v.data() == nullptr);
  }

  SECTION("SizeConstructorZero") {
    DeviceVector<int> v(0);
    CHECK(v.empty() == true);
    CHECK(v.size() == 0);
    CHECK(v.capacity() == 0);
    CHECK_NOTHROW(v.clear());
    CHECK(v.empty() == true);
    CHECK(v.data() == nullptr);
  }

  SECTION("SizeConstructorNonzero") {
    constexpr std::size_t n = 5;
    DeviceVector<int> v(n);
    CHECK(v.empty() == false);
    CHECK(v.size() == n);
    CHECK(v.capacity() == n);
    CHECK(v.data() != nullptr);

    const std::vector<int> expected(n, 4);
    apo_test::CopyToDevice<int>(v, expected);
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("HostVectorConstructor") {
    std::vector<int> expected = {1, 2, 3, 4};
    DeviceVector<int> v(expected);

    CHECK(v.empty() == false);
    CHECK(v.size() == expected.size());
    CHECK(v.capacity() == expected.size());
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("HostRvalueConstructor") {
    const std::vector<int> expected = {1, 2, 3, 4};
    DeviceVector<int> v(std::vector<int>{1, 2, 3, 4});

    CHECK(v.empty() == false);
    CHECK(v.size() == expected.size());
    CHECK(v.capacity() == expected.size());
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("CopyConstructor") {
    const std::vector<int> expected = {1, 2, 3, 4, 5};
    const std::vector<int> changed = {-1, 2, 3, 4, 5};

    DeviceVector<int> v1(expected);
    DeviceVector<int> v2(v1);

    CHECK(v1.size() == v2.size());
    CHECK(v1.capacity() == v2.capacity());
    CHECK(v1.data() != v2.data());

    apo_test::CopyToDevice<int>(v1, changed);
    CHECK(apo_test::CopyToHost<int>(v1) == changed);
    CHECK(apo_test::CopyToHost<int>(v2) == expected);
  }
}

TEST_CASE("CapacityAndResizeBehavior") {
  SECTION("ResizeFromEmpty") {
    constexpr std::size_t n = 5;
    DeviceVector<int> v;
    v.resize(n);

    CHECK(v.empty() == false);
    CHECK(v.size() == n);
    CHECK(v.capacity() == n);
    CHECK(v.data() != nullptr);
  }

  SECTION("ResizeGrowPreserve") {
    DeviceVector<int> v({1, 2, 3});
    v.resize(8);

    CHECK(v.size() == 8);
    CHECK(v.capacity() == 8);

    const std::vector<int> expectedPrefix = {1, 2, 3};
    CHECK(apo_test::CopyToHost<int>(v, expectedPrefix.size()) ==
          expectedPrefix);
  }

  SECTION("ResizeShrinkPreserve") {
    DeviceVector<int> v({1, 2, 3, 4, 5, 6});
    v.resize(3);

    const std::vector<int> expected = {1, 2, 3};
    CHECK(v.size() == 3);
    CHECK(v.capacity() == 6);
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("ResizeToZeroKeepsCapacity") {
    DeviceVector<int> v({1, 2, 3, 4});
    v.resize(0);

    CHECK(v.empty() == true);
    CHECK(v.size() == 0);
    CHECK(v.capacity() == 4);
    CHECK(v.data() != nullptr);
  }

  SECTION("ShrinkToFitNonempty") {
    DeviceVector<int> v({1, 2, 3, 4, 5, 6});
    v.resize(3);
    v.shrink_to_fit();

    const std::vector<int> expected = {1, 2, 3};
    CHECK(v.size() == 3);
    CHECK(v.capacity() == 3);
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("ShrinkToFitEmpty") {
    DeviceVector<int> v;
    CHECK_NOTHROW(v.shrink_to_fit());
    CHECK(v.empty() == true);
    CHECK(v.size() == 0);
    CHECK(v.capacity() == 0);
    CHECK(v.data() == nullptr);
  }

  SECTION("Clear") {
    DeviceVector<int> v({1, 2, 3, 4, 5, 6});
    v.clear();

    CHECK(v.empty() == true);
    CHECK(v.size() == 0);
    CHECK(v.capacity() == 0);
    CHECK(v.data() == nullptr);

    CHECK_NOTHROW(v.clear());
    CHECK(v.empty() == true);
    CHECK(v.size() == 0);
    CHECK(v.capacity() == 0);
    CHECK(v.data() == nullptr);
  }

  SECTION("PushBackFromEmpty") {
    DeviceVector<int> v;
    v.push_back(4);

    CHECK(v.size() == 1);
    CHECK(v.capacity() >= 1);

    const std::vector<int> expected = {4};
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("PushBackSecondElement") {
    DeviceVector<int> v(std::vector<int>{1});
    v.push_back(4);

    CHECK(v.size() == 2);
    CHECK(v.capacity() >= 2);

    const std::vector<int> expected = {1, 4};
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("PushBackMany") {
    std::vector<int> expected(100);
    for (std::size_t i = 0; i < expected.size(); i++)
      expected[i] = static_cast<int>(i + 1);

    DeviceVector<int> v;
    for (std::size_t i = 0; i < expected.size(); i++)
      v.push_back(static_cast<int>(i + 1));

    CHECK(v.size() == expected.size());
    CHECK(v.capacity() >= expected.size());
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("PushBackAfterClear") {
    DeviceVector<int> v({1, 2, 3, 4});
    v.clear();
    v.push_back(64);

    CHECK(v.size() == 1);
    CHECK(v.capacity() >= 1);

    const std::vector<int> expected = {64};
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("PushBackAfterShrink") {
    DeviceVector<int> v(8);
    const std::vector<int> prefix = {1, 2, 3, 4};
    apo_test::CopyToDevice<int>(v, prefix);

    v.resize(4);
    v.shrink_to_fit();
    v.push_back(64);

    const std::vector<int> expected = {1, 2, 3, 4, 64};
    CHECK(v.size() == expected.size());
    CHECK(v.capacity() >= expected.size());
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }
}

TEST_CASE("CopyAssignmentSwap") {
  SECTION("VectorAssign") {
    const std::vector<int> expected = {1, 2, 3, 4};
    DeviceVector<int> v;

    v = expected;

    CHECK(v.size() == expected.size());
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("RvalueVectorAssign") {
    const std::vector<int> expected = {5, 6, 7};
    DeviceVector<int> v;

    v = std::vector<int>{5, 6, 7};

    CHECK(v.size() == expected.size());
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("CopyAssignmentDeepCopy") {
    const std::vector<int> expected = {1, 2, 3, 4};
    const std::vector<int> changed = {4, 3, 2, 1};

    DeviceVector<int> v1(expected);
    DeviceVector<int> v2(100);
    DeviceVector<int> v3(2);

    v2 = v1;
    v3 = v1;

    CHECK(v1.size() == v2.size());
    CHECK(v1.size() == v3.size());
    CHECK(v1.capacity() == v2.capacity());
    CHECK(v1.capacity() == v3.capacity());
    CHECK(v1.data() != v2.data());
    CHECK(v1.data() != v3.data());

    apo_test::CopyToDevice<int>(v1, changed);

    CHECK(apo_test::CopyToHost<int>(v1) == changed);
    CHECK(apo_test::CopyToHost<int>(v2) == expected);
    CHECK(apo_test::CopyToHost<int>(v3) == expected);
  }

  SECTION("RvalueConstructor") {
    const std::vector<int> expected = {1, 2, 3, 4};
    DeviceVector<int> v1(expected);
    DeviceVector<int> v2(std::move(v1));

    CHECK(v2.size() == expected.size());
    CHECK(apo_test::CopyToHost<int>(v2) == expected);
  }

  SECTION("RvalueAssignment") {
    const std::vector<int> expected = {1, 2, 3, 4};
    DeviceVector<int> v1(expected);
    DeviceVector<int> v2({9});

    v2 = std::move(v1);

    CHECK(v2.size() == expected.size());
    CHECK(apo_test::CopyToHost<int>(v2) == expected);
  }

  SECTION("SelfAssignment") {
    const std::vector<int> expected = {1, 2, 3, 4};
    DeviceVector<int> v(expected);

    CHECK_NOTHROW(v = v);
    CHECK(v.size() == expected.size());
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("SelfMoveAssignment") {
    const std::vector<int> expected = {1, 2, 3, 4};
    DeviceVector<int> v(expected);

    CHECK_NOTHROW(v = std::move(v));
    CHECK(v.size() == expected.size());
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }

  SECTION("SwapNonemptyNonempty") {
    DeviceVector<int> v1({1, 2, 3, 4});
    DeviceVector<int> v2({5, 4, 3, 2, 1});

    CHECK_NOTHROW(v1.swap(v2));

    CHECK(v1.size() == 5);
    CHECK(v1.capacity() == 5);
    CHECK(v2.size() == 4);
    CHECK(v2.capacity() == 4);

    const std::vector<int> expectedV1 = {5, 4, 3, 2, 1};
    const std::vector<int> expectedV2 = {1, 2, 3, 4};
    CHECK(apo_test::CopyToHost<int>(v1) == expectedV1);
    CHECK(apo_test::CopyToHost<int>(v2) == expectedV2);
  }

  SECTION("SwapEmptyNonempty") {
    DeviceVector<int> v1;
    DeviceVector<int> v2({1, 2, 3, 4});

    CHECK_NOTHROW(v1.swap(v2));

    CHECK(v1.empty() == false);
    CHECK(v1.size() == 4);
    CHECK(v1.capacity() == 4);
    CHECK(v2.empty() == true);
    CHECK(v2.size() == 0);
    CHECK(v2.capacity() == 0);

    const std::vector<int> expected = {1, 2, 3, 4};
    CHECK(apo_test::CopyToHost<int>(v1) == expected);
  }

  SECTION("SwapNonemptyEmpty") {
    DeviceVector<int> v1({1, 2, 3, 4});
    DeviceVector<int> v2;

    CHECK_NOTHROW(v1.swap(v2));

    CHECK(v1.empty() == true);
    CHECK(v1.size() == 0);
    CHECK(v1.capacity() == 0);
    CHECK(v2.empty() == false);
    CHECK(v2.size() == 4);
    CHECK(v2.capacity() == 4);

    const std::vector<int> expected = {1, 2, 3, 4};
    CHECK(apo_test::CopyToHost<int>(v2) == expected);
  }

  SECTION("SwapSelf") {
    const std::vector<int> expected = {1, 2, 3, 4};
    DeviceVector<int> v(expected);

    CHECK_NOTHROW(v.swap(v));
    CHECK(v.size() == 4);
    CHECK(v.capacity() == 4);
    CHECK(apo_test::CopyToHost<int>(v) == expected);
  }
}
