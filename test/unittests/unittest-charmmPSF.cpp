// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "CharmmPSF.h"
#include "apo_test_helpers.h"
#include "catch.hpp"

TEST_CASE("CharmPSFDefaultConstructor") {
  CharmmPSF psf;

  CHECK(psf.getNumAtoms() == -1);
  CHECK(psf.getNumBonds() == -1);
  CHECK(psf.getNumAngles() == -1);
  CHECK(psf.getNumDihedrals() == -1);
  CHECK(psf.getNumImpropers() == -1);
  CHECK(psf.getNumCrossTerms() == -1);

  CHECK(psf.getSegmentIdentifiers().empty() == true);
  CHECK(psf.getResidueIdentifiers().empty() == true);
  CHECK(psf.getResidueNames().empty() == true);
  CHECK(psf.getAtomNames().empty() == true);
  CHECK(psf.getAtomTypes().empty() == true);
  CHECK(psf.getCharges().empty() == true);
  CHECK(psf.getMasses().empty() == true);

  CHECK(psf.getBonds().empty() == true);
  CHECK(psf.getAngles().empty() == true);
  CHECK(psf.getDihedrals().empty() == true);
  CHECK(psf.getImpropers().empty() == true);
  CHECK(psf.getCrossTerms().empty() == true);

  CHECK(psf.getConnected12().empty() == true);
  CHECK(psf.getConnected13().empty() == true);
  CHECK(psf.getConnected14().empty() == true);
  CHECK(psf.getIblo14().empty() == true);
  CHECK(psf.getInb14().empty() == true);

  CHECK(psf.getWaterMolecules().size() == 0);
  CHECK(psf.getResidues().size() == 0);
  CHECK(psf.getGroups().size() == 0);
  CHECK(psf.getFileName().empty() == true);
}

TEST_CASE("CharmmPSFParsesLinearChain") {
  const std::string fileName = "tmp_charmm_psf_linear_chain.psf";
  const std::string psfText = R"PSF(PSF

       1 !NTITLE
 REMARKS generated CharmmPSF linear-chain unit test
       4 !NATOM
       1 SEG1     1 RES1 C1   CT1   -0.300000  12.0110           0
       2 SEG1     1 RES1 H1   HA1    0.100000   1.0080           0
       3 SEG1     2 RES2 C2   CT2    0.100000  12.0110           0
       4 SEG1     2 RES2 H2   HA2    0.100000   1.0080           0
       3 !NBOND: bonds
       1       2       2       3       3       4
       2 !NTHETA: angles
       1       2       3       2       3       4
       1 !NPHI: dihedrals
       1       2       3       4
       0 !NIMPHI: impropers
       0 !NDON: donors
       0 !NACC: acceptors
       0 !NCRTERM: cross-terms
)PSF";

  apo_test::WriteTextFile(fileName, psfText);

  CharmmPSF psf(fileName);

  SECTION("CountsAndFileName") {
    CHECK(psf.getFileName() == fileName);
    CHECK(psf.getNumAtoms() == 4);
    CHECK(psf.getNumBonds() == 3);
    CHECK(psf.getNumAngles() == 2);
    CHECK(psf.getNumDihedrals() == 1);
    CHECK(psf.getNumImpropers() == 0);
    CHECK(psf.getNumCrossTerms() == 0);
  }

  SECTION("AtomMetadata") {
    CHECK(psf.getSegmentIdentifiers() ==
          std::vector<std::string>{"SEG1", "SEG1", "SEG1", "SEG1"});
    CHECK(psf.getResidueIdentifiers() == std::vector<int>{1, 1, 2, 2});
    CHECK(psf.getResidueNames() ==
          std::vector<std::string>{"RES1", "RES1", "RES2", "RES2"});
    CHECK(psf.getAtomNames() ==
          std::vector<std::string>{"C1", "H1", "C2", "H2"});
    CHECK(psf.getAtomTypes() ==
          std::vector<std::string>{"CT1", "HA1", "CT2", "HA2"});

    REQUIRE(psf.getCharges().size() == 4);
    CHECK(psf.getCharges()[0] == Approx(-0.3));
    CHECK(psf.getCharges()[1] == Approx(0.1));
    CHECK(psf.getCharges()[2] == Approx(0.1));
    CHECK(psf.getCharges()[3] == Approx(0.1));

    REQUIRE(psf.getMasses().size() == 4);
    CHECK(psf.getMasses()[0] == Approx(12.011));
    CHECK(psf.getMasses()[1] == Approx(1.008));
    CHECK(psf.getMasses()[2] == Approx(12.011));
    CHECK(psf.getMasses()[3] == Approx(1.008));

    CHECK(psf.getNetCharge() == Approx(0.0).margin(1.0e-12));
    CHECK(psf.getTotalMass() == Approx(26.038));
  }

  SECTION("BondedTerms") {
    REQUIRE(psf.getBonds().size() == 3);
    CHECK(psf.getBonds()[0].iatom == 0);
    CHECK(psf.getBonds()[0].jatom == 1);
    CHECK(psf.getBonds()[1].iatom == 1);
    CHECK(psf.getBonds()[1].jatom == 2);
    CHECK(psf.getBonds()[2].iatom == 2);
    CHECK(psf.getBonds()[2].jatom == 3);

    REQUIRE(psf.getAngles().size() == 2);
    CHECK(psf.getAngles()[0].iatom == 0);
    CHECK(psf.getAngles()[0].jatom == 1);
    CHECK(psf.getAngles()[0].katom == 2);
    CHECK(psf.getAngles()[1].iatom == 1);
    CHECK(psf.getAngles()[1].jatom == 2);
    CHECK(psf.getAngles()[1].katom == 3);

    REQUIRE(psf.getDihedrals().size() == 1);
    CHECK(psf.getDihedrals()[0].iatom == 0);
    CHECK(psf.getDihedrals()[0].jatom == 1);
    CHECK(psf.getDihedrals()[0].katom == 2);
    CHECK(psf.getDihedrals()[0].latom == 3);

    CHECK(psf.getImpropers().empty() == true);
    CHECK(psf.getCrossTerms().empty() == true);
  }

  SECTION("ResiduesGroupsAndWater") {
    const std::vector<int2> &residues = psf.getResidues().getHostArray();
    REQUIRE(residues.size() == 2);
    CHECK(residues[0].x == 0);
    CHECK(residues[0].y == 1);
    CHECK(residues[1].x == 2);
    CHECK(residues[1].y == 3);

    const std::vector<int2> &groups = psf.getGroups().getHostArray();
    REQUIRE(groups.size() == 1);
    CHECK(groups[0].x == 0);
    CHECK(groups[0].y == 3);

    CHECK(psf.getWaterMolecules().size() == 0);
    CHECK(psf.getWaterMolecules().getHostArray().empty() == true);
  }

  SECTION("ConnectedComponentsAndExclusions") {
    const std::vector<std::set<int>> expected12 = {
        std::set<int>{1}, std::set<int>{0, 2}, std::set<int>{1, 3},
        std::set<int>{2}};
    const std::vector<std::set<int>> expected13 = {
        std::set<int>{2}, std::set<int>{3}, std::set<int>{0}, std::set<int>{1}};
    const std::vector<std::set<int>> expected14 = {
        std::set<int>{3}, std::set<int>{}, std::set<int>{}, std::set<int>{0}};

    CHECK(psf.getConnected12() == expected12);
    CHECK(psf.getConnected13() == expected13);
    CHECK(psf.getConnected14() == expected14);

    CHECK(psf.getIblo14() == std::vector<int>{3, 5, 6, 6});
    CHECK(psf.getInb14() == std::vector<int>{2, 3, 4, 3, 4, 4});

    const InclusionExclusion exclusionLists = psf.getInclusionExclusionLists();

    CHECK(exclusionLists.sizes == std::vector<int>{1, 5});
    CHECK(exclusionLists.in14_ex14 ==
          std::vector<int>{0, 3, 0, 1, 0, 2, 1, 2, 1, 3, 2, 3});
  }

  SECTION("CopyConstructorDeepCopy") {
    CharmmPSF copy(psf);

    CHECK(copy.getFileName() == psf.getFileName());
    CHECK(copy.getNumAtoms() == psf.getNumAtoms());
    CHECK(copy.getNumBonds() == psf.getNumBonds());
    CHECK(copy.getCharges().data() != psf.getCharges().data());
    CHECK(copy.getAtomNames().data() != psf.getAtomNames().data());

    copy.getCharges()[0] = 9.0;
    copy.getAtomNames()[0] = "XX";

    CHECK(psf.getCharges()[0] == Approx(-0.3));
    CHECK(psf.getAtomNames()[0] == "C1");
    CHECK(copy.getCharges()[0] == Approx(9.0));
    CHECK(copy.getAtomNames()[0] == "XX");
    CHECK(copy.getNetCharge() == Approx(9.3));
  }

  SECTION("RvalueConstructor") {
    CharmmPSF source(fileName);
    CharmmPSF copy(std::move(source));

    CHECK(copy.getFileName() == fileName);
    CHECK(copy.getNumAtoms() == 4);
    CHECK(copy.getNumBonds() == 3);
    CHECK(copy.getAtomNames() ==
          std::vector<std::string>{"C1", "H1", "C2", "H2"});
    CHECK(copy.getNetCharge() == Approx(0.0).margin(1.0e-12));
  }

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("CharmmPSFDetectsWaterMolecules") {
  const std::string fileName = "tmp_charmm_psf_water.psf";
  const std::string psfText = R"PSF(PSF

       1 !NTITLE
 REMARKS generated CharmmPSF water unit test
       3 !NATOM
       1 WAT      1 TIP3 OH2  OT    -0.834000  15.9994           0
       2 WAT      1 TIP3 H1   HT     0.417000   1.0080           0
       3 WAT      1 TIP3 H2   HT     0.417000   1.0080           0
       2 !NBOND: bonds
       1       2       1       3
       1 !NTHETA: angles
       2       1       3
       0 !NPHI: dihedrals
       0 !NIMPHI: impropers
       0 !NDON: donors
       0 !NACC: acceptors
       0 !NCRTERM: cross-terms
)PSF";

  apo_test::WriteTextFile(fileName, psfText);

  CharmmPSF psf(fileName);

  CHECK(psf.getNumAtoms() == 3);
  CHECK(psf.getNumBonds() == 2);
  CHECK(psf.getNumAngles() == 1);
  CHECK(psf.getNetCharge() == Approx(0.0).margin(1.0e-12));

  const std::vector<int4> &waterMolecules =
      psf.getWaterMolecules().getHostArray();
  REQUIRE(waterMolecules.size() == 1);
  CHECK(waterMolecules[0].x == 0);
  CHECK(waterMolecules[0].y == 1);
  CHECK(waterMolecules[0].z == 2);
  CHECK(waterMolecules[0].w == 0);

  const std::vector<int2> &residues = psf.getResidues().getHostArray();
  REQUIRE(residues.size() == 1);
  CHECK(residues[0].x == 0);
  CHECK(residues[0].y == 2);

  const std::vector<int2> &groups = psf.getGroups().getHostArray();
  REQUIRE(groups.size() == 1);
  CHECK(groups[0].x == 0);
  CHECK(groups[0].y == 2);

  apo_test::RemoveIfExists(fileName);
}

TEST_CASE("CharmmPSFMalformedInpuThrows") {
  const std::string fileName = "tmp_charmm_psf_malformed.psf";
  const std::string psfText = R"PSF(PSF

       1 !NTITLE
 REMARKS missing required sections
       1 !NATOM
       1 SEG1     1 RES1 C1   CT1    0.000000  12.0110           0
)PSF";

  apo_test::WriteTextFile(fileName, psfText);

  CHECK_THROWS_AS(CharmmPSF(fileName), std::runtime_error);

  apo_test::RemoveIfExists(fileName);
}
