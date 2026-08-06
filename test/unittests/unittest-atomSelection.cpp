// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "ApoCharmmError.h"
#include "AtomSelection.h"
#include "AtomSelector.h"
#include "CharmmPSF.h"
#include "apo_test_helpers.h"
#include "catch.hpp"

#include <limits>
#include <string>
#include <string_view>

namespace {

const std::string TEST_PSF_TEXT = R"PSF(PSF

       1 !NTITLE
 REMARKS generated AtomSelector unit test PSF
      10 !NATOM
       1 SEG1     1 ALA  N    NH1   -0.300000  14.0070           0
       2 SEG1     1 ALA  CA   CT1    0.100000  12.0110           0
       3 SEG1     1 ALA  CB   CT2    0.000000  12.0110           0
       4 SEG1     1 ALA  HB1  HA     0.100000   1.0080           0
       5 SEG2     2 GLY  N    NH1   -0.200000  14.0070           0
       6 SEG2     2 GLY  CA   CT1    0.100000  12.0110           0
       7 SEG2     2 GLY  HA1  HA     0.100000   1.0080           0
       8 WAT      3 TIP3 OH2  OT    -0.834000  15.9994           0
       9 WAT      3 TIP3 H1   HT     0.417000   1.0080           0
      10 WAT      3 TIP3 H2   HT     0.417000   1.0080           0
       7 !NBOND: bonds
       1       2       2       3       3       4       5       6
       6       7       8       9       8      10
       4 !NTHETA: angles
       1       2       3       2       3       4       5       6       7
       9       8      10
       1 !NPHI: dihedrals
       1       2       3       4
       0 !NIMPHI: impropers
       0 !NDON: donors
       0 !NACC: acceptors
       0 !NCRTERM: cross-terms
)PSF";

void CheckSelection(const AtomSelection &selection, const int numAtoms,
                    const std::vector<int> &expectedIndices) {
  CHECK(selection.getNumAtoms() == numAtoms);
  CHECK(selection.getNumSelected() == static_cast<int>(expectedIndices.size()));
  CHECK(selection.getAtomIndices() == expectedIndices);

  std::vector<bool> expectedContains(static_cast<std::size_t>(numAtoms), false);
  for (const int index : expectedIndices)
    expectedContains[index] = true;

  for (int i = 0; i < numAtoms; i++)
    CHECK(selection.contains(i) == expectedContains[i]);

  return;
}

} // namespace

TEST_CASE("AtomSelectionBitsetBehavior") {
  SECTION("ConstructNone") {
    AtomSelection selection(130);

    CHECK(selection.getNumAtoms() == 130);
    CHECK(selection.getNumSelected() == 0);
    CHECK(selection.getAtomIndices().empty() == true);
    CHECK(selection.contains(0) == false);
    CHECK(selection.contains(63) == false);
    CHECK(selection.contains(64) == false);
    CHECK(selection.contains(129) == false);
  }

  SECTION("ConstructAllMasksUnusedBits") {
    AtomSelection selection(130, AtomSelection::InitialValue::ALL);

    CHECK(selection.getNumAtoms() == 130);
    CHECK(selection.getNumSelected() == 130);
    CHECK(selection.getAtomIndices() == apo_test::MakeRange(0, 130));
    CHECK(selection.contains(0) == true);
    CHECK(selection.contains(63) == true);
    CHECK(selection.contains(64) == true);
    CHECK(selection.contains(129) == true);
  }

  SECTION("SetClearAndFillAcrossWordBoundaries") {
    AtomSelection selection(130);

    selection.set(0);
    selection.set(63);
    selection.set(64);
    selection.set(129);

    CheckSelection(selection, 130, {0, 63, 64, 129});

    selection.set(64, false);
    CheckSelection(selection, 130, {0, 63, 129});

    selection.clear();
    CheckSelection(selection, 130, {});

    selection.fill();
    CheckSelection(selection, 130, apo_test::MakeRange(0, 130));
  }

  SECTION("SetNumAtomsResetsSelection") {
    AtomSelection selection(10, AtomSelection::InitialValue::ALL);

    selection.setNumAtoms(65, AtomSelection::InitialValue::NONE);
    CheckSelection(selection, 65, {});

    selection.setNumAtoms(65, AtomSelection::InitialValue::ALL);
    CheckSelection(selection, 65, apo_test::MakeRange(0, 65));
  }

  SECTION("CopyConstructorAndAssignmentDeepCopy") {
    AtomSelection source(70);
    source.set(0);
    source.set(64);
    source.set(69);

    AtomSelection copy(source);
    AtomSelection assigned(70);
    assigned = source;

    source.set(64, false);

    CheckSelection(source, 70, {0, 69});
    CheckSelection(copy, 70, {0, 64, 69});
    CheckSelection(assigned, 70, {0, 64, 69});
  }

  SECTION("RvalueConstructorAndAssignment") {
    AtomSelection source(70);
    source.set(1);
    source.set(65);

    AtomSelection moved(std::move(source));
    CheckSelection(moved, 70, {1, 65});

    AtomSelection assigned(70);
    assigned = std::move(moved);
    CheckSelection(assigned, 70, {1, 65});
  }

  SECTION("LogicalOperators") {
    AtomSelection left(130);
    left.set(0);
    left.set(2);
    left.set(64);

    AtomSelection right(130);
    right.set(2);
    right.set(3);
    right.set(64);
    right.set(129);

    AtomSelection intersection = left & right;
    AtomSelection unionSelection = left | right;

    CheckSelection(intersection, 130, {2, 64});
    CheckSelection(unionSelection, 130, {0, 2, 3, 64, 129});

    left &= right;
    CheckSelection(left, 130, {2, 64});

    right |= intersection;
    CheckSelection(right, 130, {2, 3, 64, 129});
  }

  SECTION("RejectsInvalidConstructionAndIndices") {
    apo_test::CheckApoCharmmError(
        []() -> void { (void)AtomSelection(-1); },
        ApoCharmmErrorCode::InvalidArgument,
        "Atom count must be non-negative; observed -1");

    const int minimumAtomCount = std::numeric_limits<int>::min();
    apo_test::CheckApoCharmmError(
        [minimumAtomCount]() -> void { (void)AtomSelection(minimumAtomCount); },
        ApoCharmmErrorCode::InvalidArgument,
        "Atom count must be non-negative; observed " +
            std::to_string(minimumAtomCount));

    AtomSelection selection(4);
    selection.set(2);

    apo_test::CheckApoCharmmError(
        [&selection]() -> void {
          selection.setNumAtoms(-1, AtomSelection::InitialValue::ALL);
        },
        ApoCharmmErrorCode::InvalidArgument,
        "Atom count must be non-negative; observed -1");
    CheckSelection(selection, 4, {2});

    apo_test::CheckApoCharmmError(
        [&selection]() -> void { (void)selection.contains(-1); },
        ApoCharmmErrorCode::InvalidArgument,
        "Atom index is out of range; expected [0, 4), observed -1");
    apo_test::CheckApoCharmmError(
        [&selection]() -> void { (void)selection.contains(4); },
        ApoCharmmErrorCode::InvalidArgument,
        "Atom index is out of range; expected [0, 4), observed 4");
    apo_test::CheckApoCharmmError(
        [&selection]() -> void { selection.set(-1); },
        ApoCharmmErrorCode::InvalidArgument,
        "Atom index is out of range; expected [0, 4), observed -1");
    apo_test::CheckApoCharmmError(
        [&selection]() -> void { selection.set(4); },
        ApoCharmmErrorCode::InvalidArgument,
        "Atom index is out of range; expected [0, 4), observed 4");
    CheckSelection(selection, 4, {2});
  }

  SECTION("RejectsIncompatibleLogicalOperations") {
    AtomSelection left(4);
    AtomSelection right(5);

    apo_test::CheckApoCharmmError(
        [&left, &right]() -> void { left &= right; },
        ApoCharmmErrorCode::InvalidArgument,
        "Selection atom count mismatch; expected 4, observed 5");
    apo_test::CheckApoCharmmError(
        [&left, &right]() -> void { left |= right; },
        ApoCharmmErrorCode::InvalidArgument,
        "Selection atom count mismatch; expected 4, observed 5");
    apo_test::CheckApoCharmmError(
        [&left, &right]() -> void { (void)(left & right); },
        ApoCharmmErrorCode::InvalidArgument,
        "Selection atom count mismatch; expected 4, observed 5");
    apo_test::CheckApoCharmmError(
        [&left, &right]() -> void { (void)(left | right); },
        ApoCharmmErrorCode::InvalidArgument,
        "Selection atom count mismatch; expected 4, observed 5");
  }
}

TEST_CASE("AtomSelectorSelectionLanguage") {
  const std::string fileName = "tmp_atom_selection_test.psf";
  apo_test::WriteTextFile(fileName, TEST_PSF_TEXT);

  auto psf = std::make_shared<CharmmPSF>(fileName);
  AtomSelector selector(psf);

  SECTION("BasicSelections") {
    CheckSelection(selector.select("none"), 10, {});
    CheckSelection(selector.select("all"), 10, apo_test::MakeRange(0, 10));

    CheckSelection(selector.select("bynum 2:4"), 10, {1, 2, 3});
    CheckSelection(selector.select("bynu 4:2"), 10, {1, 2, 3});
    CheckSelection(selector.select("bynum -5:2"), 10, {0, 1});
    CheckSelection(selector.select("bynum 20:25"), 10, {});

    CheckSelection(selector.select("type CA"), 10, {1, 5});
    CheckSelection(selector.select("type ca"), 10, {1, 5});
    CheckSelection(selector.select("chemical CT1"), 10, {1, 5});
    CheckSelection(selector.select("chemical CT%"), 10, {1, 2, 5});

    CheckSelection(selector.select("segid SEG1"), 10, {0, 1, 2, 3});
    CheckSelection(selector.select("segi SEG2"), 10, {4, 5, 6});
    CheckSelection(selector.select("resid 1"), 10, {0, 1, 2, 3});
    CheckSelection(selector.select("resi 1:2"), 10, {0, 1, 2, 3, 4, 5, 6});
    CheckSelection(selector.select("resname GLY"), 10, {4, 5, 6});
    CheckSelection(selector.select("resn TIP3"), 10, {7, 8, 9});
    CheckSelection(selector.select("atom SEG1 1 CA"), 10, {1});
  }

  SECTION("Wildcards") {
    CheckSelection(selector.select("type H*"), 10, {3, 6, 8, 9});
    CheckSelection(selector.select("type H+"), 10, {8, 9});
    CheckSelection(selector.select("type H#"), 10, {8, 9});
    CheckSelection(selector.select("type H%"), 10, {8, 9});
    CheckSelection(selector.select("type C%"), 10, {1, 2, 5});
    CheckSelection(selector.select("type C*"), 10, {1, 2, 5});
  }

  SECTION("LogicalOperatorsAndPrecedence") {
    CheckSelection(selector.select("type CA .and. resn ALA"), 10, {1});
    CheckSelection(selector.select("type CA .or. type N .and. resn ALA"), 10,
                   {0, 1, 5});
    CheckSelection(selector.select("(type CA .or. type N) .and. resn ALA"), 10,
                   {0, 1});
    CheckSelection(selector.select(".not. type H*"), 10, {0, 1, 2, 4, 5, 7});
  }

  SECTION("PrefixExpansionOperators") {
    CheckSelection(selector.select(".byres. type HB1"), 10, {0, 1, 2, 3});
    CheckSelection(selector.select(".byres. (type CA .and. resn GLY)"), 10,
                   {4, 5, 6});
    CheckSelection(selector.select(".bygroup. type HA1"), 10, {4, 5, 6});
    CheckSelection(selector.select(".bygroup. type H1"), 10, {7, 8, 9});

    CheckSelection(selector.select(".bonded. (atom SEG1 1 CA)"), 10, {0, 2});
    CheckSelection(selector.select(".bonded. bynum 8"), 10, {8, 9});
  }

  SECTION("InvalidSelectionsThrow") {
    CHECK_THROWS(selector.select(""));
    CHECK_THROWS(selector.select(".and. type CA"));
    CHECK_THROWS(selector.select("type CA resn ALA"));
    CHECK_THROWS(selector.select("(type CA .or. type N"));
    CHECK_THROWS(selector.select("type CA .or."));
    CHECK_THROWS(selector.select("bynu A:C"));
    CHECK_THROWS(selector.select(".around. type CA"));
  }

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("AtomSelectorRejectsNullPSF") {
  AtomSelector selector(nullptr);

  CHECK_THROWS_AS(selector.select("all"), std::runtime_error);
}
