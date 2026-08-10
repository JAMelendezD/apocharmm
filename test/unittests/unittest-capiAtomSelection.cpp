// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "AtomSelection.h"
#include "apo_test_helpers.h"
#include "apocharmm_c/AtomSelection.h"
#include "apocharmm_c/Error.h"
#include "apocharmm_c/Status.h"
#include "apocharmm_c/detail/AtomSelectionHandle.h"
#include "catch.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct AtomSelectionDeleter {
  void operator()(apo_atom_selection *selection) const noexcept {
    apo_atom_selection_destroy(selection);
    return;
  }
};

using AtomSelectionHandle =
    std::unique_ptr<apo_atom_selection, AtomSelectionDeleter>;

AtomSelectionHandle MakeSelection(const int numAtoms,
                                  const std::vector<int> &selectedAtoms) {
  AtomSelectionHandle handle(new apo_atom_selection());
  handle->object = std::make_shared<AtomSelection>(numAtoms);

  for (const int atomIndex : selectedAtoms)
    handle->object->set(atomIndex);

  return handle;
}

} // namespace

TEST_CASE("CapiAtomSelectionAccessors") {
  AtomSelectionHandle selection = MakeSelection(5, {1, 3});

  std::size_t numAtoms = 0;
  CHECK(apo_atom_selection_get_num_atoms(&numAtoms, selection.get()) ==
        APO_STATUS_OK);
  CHECK(numAtoms == 5);

  std::size_t numSelected = 0;
  CHECK(apo_atom_selection_get_num_selected(&numSelected, selection.get()) ==
        APO_STATUS_OK);
  CHECK(numSelected == 2);

  int atomIndices[2] = {-1, -1};
  CHECK(apo_atom_selection_get_atom_indices(atomIndices, 2, selection.get()) ==
        APO_STATUS_OK);
  CHECK(atomIndices[0] == 1);
  CHECK(atomIndices[1] == 3);

  bool isSelected = false;
  CHECK(apo_atom_selection_contains(&isSelected, selection.get(), 1) ==
        APO_STATUS_OK);
  CHECK(isSelected == true);

  CHECK(apo_atom_selection_contains(&isSelected, selection.get(), 2) ==
        APO_STATUS_OK);
  CHECK(isSelected == false);

  AtomSelectionHandle emptySelection = MakeSelection(0, {});
  CHECK(apo_atom_selection_get_atom_indices(nullptr, 0, emptySelection.get()) ==
        APO_STATUS_OK);

  CHECK(std::string(apo_last_error()).empty() == true);
}

TEST_CASE("CapiAtomSelectionValidatesPointersAndBuffers") {
  SECTION("GetNumAtoms") {
    AtomSelectionHandle selection = MakeSelection(4, {1, 3});

    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_get_num_atoms(nullptr, selection.get()),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_get_num_atoms: num_atoms is NULL");

    std::size_t numAtoms = 17;
    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_get_num_atoms(&numAtoms, nullptr),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_get_num_atoms: AtomSelection is NULL");
    CHECK(numAtoms == 0);

    apo_atom_selection emptyHandle;
    numAtoms = 17;
    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_get_num_atoms(&numAtoms, &emptyHandle),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_get_num_atoms: AtomSelection object is NULL");
    CHECK(numAtoms == 0);
  }

  SECTION("GetNumSelected") {
    AtomSelectionHandle selection = MakeSelection(4, {1, 3});

    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_get_num_selected(nullptr, selection.get()),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_get_num_selected: num_selected is NULL");

    std::size_t numSelected = 17;
    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_get_num_selected(&numSelected, nullptr),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_get_num_selected: AtomSelection is NULL");
    CHECK(numSelected == 0);

    apo_atom_selection emptyHandle;
    numSelected = 17;
    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_get_num_selected(&numSelected, &emptyHandle),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_get_num_selected: AtomSelection object is NULL");
    CHECK(numSelected == 0);
  }

  SECTION("GetAtomIndices") {
    AtomSelectionHandle selection = MakeSelection(4, {1, 3});
    int atomIndices[2] = {-1, -1};

    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_get_atom_indices(atomIndices, 2, nullptr),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_get_atom_indices: AtomSelection is NULL");

    apo_atom_selection emptyHandle;
    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_get_atom_indices(atomIndices, 2, &emptyHandle),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_get_atom_indices: AtomSelection object is NULL");

    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_get_atom_indices(nullptr, 2, selection.get()),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_get_atom_indices: Atom index buffer is NULL");

    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_get_atom_indices(atomIndices, 1, selection.get()),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_get_atom_indices: Atom index buffer is too small");
  }

  SECTION("Contains") {
    AtomSelectionHandle selection = MakeSelection(4, {1, 3});

    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_contains(nullptr, selection.get(), 1),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_contains: is_selected is NULL");

    bool isSelected = true;
    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_contains(&isSelected, nullptr, 1),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_contains: AtomSelection is NULL");
    CHECK(isSelected == false);

    apo_atom_selection emptyHandle;
    isSelected = true;
    apo_test::CheckStatusAndDiagnostic(
        apo_atom_selection_contains(&isSelected, &emptyHandle, 1),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_atom_selection_contains: AtomSelection object is NULL");
    CHECK(isSelected == false);

    apo_status status = APO_STATUS_OK;

    isSelected = true;
    CHECK_NOTHROW((status = apo_atom_selection_contains(&isSelected,
                                                        selection.get(), -1)));
    apo_test::CheckNativeError(
        status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
        "apo_atom_selection_contains",
        "Atom index is out of range; expected [0, 4), observed -1",
        "src/AtomSelection.cpp", "checkIndex");
    CHECK(isSelected == false);

    isSelected = true;
    CHECK_NOTHROW((
        status = apo_atom_selection_contains(&isSelected, selection.get(), 4)));
    apo_test::CheckNativeError(
        status, APO_STATUS_INVALID_ARGUMENT, "InvalidArgument",
        "apo_atom_selection_contains",
        "Atom index is out of range; expected [0, 4), observed 4",
        "src/AtomSelection.cpp", "checkIndex");
    CHECK(isSelected == false);
  }
}

TEST_CASE("CapiAtomSelectionSuccessClearsStaleDiagnostic") {
  AtomSelectionHandle selection = MakeSelection(4, {1, 3});

  std::size_t numAtoms = 17;
  REQUIRE(apo_atom_selection_get_num_atoms(&numAtoms, nullptr) ==
          APO_STATUS_INVALID_ARGUMENT);
  REQUIRE(std::string(apo_last_error()).empty() == false);
  REQUIRE(numAtoms == 0);

  CHECK(apo_atom_selection_get_num_atoms(&numAtoms, selection.get()) ==
        APO_STATUS_OK);
  CHECK(numAtoms == 4);
  CHECK(std::string(apo_last_error()).empty() == true);
}

TEST_CASE("CapiAtomSelectionDestroyIsNoexcept") {
  AtomSelectionHandle selection = MakeSelection(4, {1, 3});
  apo_atom_selection *const rawSelection = selection.release();

  CHECK_NOTHROW(apo_atom_selection_destroy(rawSelection));
  CHECK(std::string(apo_last_error()).empty() == true);

  CHECK_NOTHROW(apo_atom_selection_destroy(nullptr));
  CHECK(std::string(apo_last_error()).empty() == true);
}
