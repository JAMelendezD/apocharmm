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
#include "CharmmParameters.h"
#include "apo_test_helpers.h"
#include "catch.hpp"

namespace {

constexpr double DOUBLE_TOLERANCE = 1.0e-12;
constexpr double FLOAT_TOLERANCE = 1.0e-6;

const std::string TEST_PARAMETER_TEXT =
    R"PRM(* generated CharmmParameters unit test file
*

BONDS
A      B      100.0     1.1
B      C      200.0     1.2
C      D      300.0     1.3

ANGLES
A      B      C       50.0    120.0    10.0    2.5
B      C      D       60.0    110.0

DIHEDRALS
A      B      C      D        0.5      3       180.0
X      B      C      X        0.2      2         0.0

IMPROPER
A      B      C      D        1.5      0       180.0

NONBONDED
A      0.0    -0.1     1.0
B      0.0    -0.2     1.1
C      0.0    -0.3     1.2
D      0.0    -0.4     1.3
E      0.0    -0.5     1.4      0.0    -0.05    1.0

NBFIX
A      D      -0.05    2.5

END
)PRM";

const std::string TEST_PARAMETER_SUPPLEMENT_TEXT =
    R"PRM(* generated supplemental CharmmParameters unit test file
*

BONDS
E      F      400.0     1.4

NONBONDED
F      0.0    -0.6      1.5

END
)PRM";

const std::string MISSING_BOND_PARAMETER_TEXT =
    R"PRM(* generated incomplete CharmmParameters unit test file
*

BONDS
A      B      100.0     1.1
B      C      200.0     1.2

END
)PRM";

const std::string MALFORMED_PARAMETER_TEXT =
    R"PRM(* generated malformed CharmmParameters unit test file
*

NONBONDED
A      0.0    -0.1     1.0      2.0

END
)PRM";

const std::string TEST_PSF_TEXT = R"PSF(PSF

       1 !NTITLE
 REMARKS generated CharmmParameters unit test PSF
       4 !NATOM
       1 SEG1     1 RES1 A1   A      0.000000  12.0110           0
       2 SEG1     1 RES1 B1   B      0.000000  12.0110           0
       3 SEG1     1 RES1 C1   C      0.000000  12.0110           0
       4 SEG1     1 RES1 D1   D      0.000000  12.0110           0
       3 !NBOND: bonds
       1       2       2       3       3       4
       2 !NTHETA: angles
       1       2       3       2       3       4
       1 !NPHI: dihedrals
       1       2       3       4
       1 !NIMPHI: impropers
       1       2       3       4
       0 !NDON: donors
       0 !NACC: acceptors
       0 !NCRTERM: cross-terms
)PSF";

} // namespace

TEST_CASE("CharmmParametersParsesSingleFile") {
  const std::string parameterFile = "tmp_charmm_parameters_full.prm";
  apo_test::WriteTextFile(parameterFile, TEST_PARAMETER_TEXT);

  CharmmParameters parameters(parameterFile);

  SECTION("OriginalFileName") {
    CHECK(parameters.getOriginalPrmFileNames() ==
          std::vector<std::string>{parameterFile});
  }

  SECTION("BondParameters") {
    const std::map<BondKey, BondValues> bonds = parameters.getBonds();

    REQUIRE(bonds.size() == 3);
    REQUIRE(bonds.count(BondKey("A", "B")) == 1);
    REQUIRE(bonds.count(BondKey("B", "C")) == 1);
    REQUIRE(bonds.count(BondKey("C", "D")) == 1);

    CHECK(bonds.at(BondKey("A", "B")).kb == Approx(100.0f));
    CHECK(bonds.at(BondKey("A", "B")).b0 == Approx(1.1f));

    CHECK(bonds.at(BondKey("B", "C")).kb == Approx(200.0f));
    CHECK(bonds.at(BondKey("B", "C")).b0 == Approx(1.2f));

    CHECK(bonds.at(BondKey("C", "D")).kb == Approx(300.0f));
    CHECK(bonds.at(BondKey("C", "D")).b0 == Approx(1.3f));
  }

  SECTION("AngleAndUreyBradleyParameters") {
    const std::map<AngleKey, AngleValues> angles = parameters.getAngles();
    const std::map<AngleKey, BondValues> ureyBradleys =
        parameters.getUreyBradleys();

    REQUIRE(angles.size() == 2);
    REQUIRE(ureyBradleys.size() == 2);

    REQUIRE(angles.count(AngleKey("A", "B", "C")) == 1);
    CHECK(angles.at(AngleKey("A", "B", "C")).kTheta == Approx(50.0f));
    CHECK(angles.at(AngleKey("A", "B", "C")).theta0 ==
          Approx(static_cast<float>(apo_test::Deg2Rad(120.0)))
              .margin(FLOAT_TOLERANCE));

    REQUIRE(angles.count(AngleKey("B", "C", "D")) == 1);
    CHECK(angles.at(AngleKey("B", "C", "D")).kTheta == Approx(60.0f));
    CHECK(angles.at(AngleKey("B", "C", "D")).theta0 ==
          Approx(static_cast<float>(apo_test::Deg2Rad(110.0)))
              .margin(FLOAT_TOLERANCE));

    CHECK(ureyBradleys.at(AngleKey("A", "B", "C")).kb == Approx(10.0f));
    CHECK(ureyBradleys.at(AngleKey("A", "B", "C")).b0 == Approx(2.5f));

    CHECK(ureyBradleys.at(AngleKey("B", "C", "D")).kb == Approx(0.0f));
    CHECK(ureyBradleys.at(AngleKey("B", "C", "D")).b0 == Approx(0.0f));
  }

  SECTION("DihedralAndImproperParameters") {
    const std::map<DihedralKey, std::vector<DihedralValues>> dihedrals =
        parameters.getDihedrals();
    const std::map<DihedralKey, ImDihedralValues> impropers =
        parameters.getImpropers();

    REQUIRE(dihedrals.count(DihedralKey("A", "B", "C", "D")) == 1);
    REQUIRE(dihedrals.at(DihedralKey("A", "B", "C", "D")).size() == 1);
    CHECK(dihedrals.at(DihedralKey("A", "B", "C", "D"))[0].kChi ==
          Approx(0.5f));
    CHECK(dihedrals.at(DihedralKey("A", "B", "C", "D"))[0].n == 3);
    CHECK(dihedrals.at(DihedralKey("A", "B", "C", "D"))[0].delta ==
          Approx(180.0f));

    REQUIRE(dihedrals.count(DihedralKey("X", "B", "C", "X")) == 1);
    REQUIRE(dihedrals.at(DihedralKey("X", "B", "C", "X")).size() == 1);
    CHECK(dihedrals.at(DihedralKey("X", "B", "C", "X"))[0].kChi ==
          Approx(0.2f));
    CHECK(dihedrals.at(DihedralKey("X", "B", "C", "X"))[0].n == 2);
    CHECK(dihedrals.at(DihedralKey("X", "B", "C", "X"))[0].delta ==
          Approx(0.0f));

    REQUIRE(impropers.count(DihedralKey("A", "B", "C", "D")) == 1);
    CHECK(impropers.at(DihedralKey("A", "B", "C", "D")).kpsi == Approx(1.5f));
    CHECK(impropers.at(DihedralKey("A", "B", "C", "D")).psi0 == Approx(180.0f));
  }

  SECTION("VdwParameters") {
    const std::map<std::string, VdwParameters> vdw =
        parameters.getVdwParameters();
    const std::map<std::string, VdwParameters> vdw14 =
        parameters.getVdw14Parameters();

    REQUIRE(vdw.size() == 5);
    REQUIRE(vdw.count("A") == 1);
    REQUIRE(vdw.count("B") == 1);
    REQUIRE(vdw.count("C") == 1);
    REQUIRE(vdw.count("D") == 1);
    REQUIRE(vdw.count("E") == 1);

    CHECK(vdw.at("A").epsilon == Approx(-0.1).margin(DOUBLE_TOLERANCE));
    CHECK(vdw.at("A").rmin_2 == Approx(1.0).margin(DOUBLE_TOLERANCE));
    CHECK(vdw.at("D").epsilon == Approx(-0.4).margin(DOUBLE_TOLERANCE));
    CHECK(vdw.at("D").rmin_2 == Approx(1.3).margin(DOUBLE_TOLERANCE));

    REQUIRE(vdw14.size() == 1);
    REQUIRE(vdw14.count("E") == 1);
    CHECK(vdw14.at("E").epsilon == Approx(-0.05).margin(DOUBLE_TOLERANCE));
    CHECK(vdw14.at("E").rmin_2 == Approx(1.0).margin(DOUBLE_TOLERANCE));
  }

  apo_test::RemoveIfExists(parameterFile);
}

TEST_CASE("CharmmParametersParsesMultipleFiles") {
  const std::string parameterFile = "tmp_charmm_parameters_full.prm";
  const std::string supplementalFile = "tmp_charmm_parameters_supplement.prm";
  apo_test::WriteTextFile(parameterFile, TEST_PARAMETER_TEXT);
  apo_test::WriteTextFile(supplementalFile,
                                TEST_PARAMETER_SUPPLEMENT_TEXT);

  CharmmParameters parameters({parameterFile, supplementalFile});

  CHECK(parameters.getOriginalPrmFileNames() ==
        std::vector<std::string>{parameterFile, supplementalFile});

  const std::map<BondKey, BondValues> bonds = parameters.getBonds();
  REQUIRE(bonds.count(BondKey("A", "B")) == 1);
  REQUIRE(bonds.count(BondKey("E", "F")) == 1);

  CHECK(bonds.at(BondKey("E", "F")).kb == Approx(400.0f));
  CHECK(bonds.at(BondKey("E", "F")).b0 == Approx(1.4f));

  const std::map<std::string, VdwParameters> vdw =
      parameters.getVdwParameters();
  REQUIRE(vdw.count("F") == 1);
  CHECK(vdw.at("F").epsilon == Approx(-0.6).margin(DOUBLE_TOLERANCE));
  CHECK(vdw.at("F").rmin_2 == Approx(1.5).margin(DOUBLE_TOLERANCE));

  apo_test::RemoveIfExists(parameterFile);
  apo_test::RemoveIfExists(supplementalFile);
}

TEST_CASE("CharmmParametersBuildsBondedParamsAndLists") {
  const std::string parameterFile = "tmp_charmm_parameters_full.prm";
  const std::string psfFile = "tmp_charmm_parameters_test.psf";
  apo_test::WriteTextFile(parameterFile, TEST_PARAMETER_TEXT);
  apo_test::WriteTextFile(psfFile, TEST_PSF_TEXT);

  CharmmParameters parameters(parameterFile);
  std::shared_ptr<CharmmPSF> psf = std::make_shared<CharmmPSF>(psfFile);

  const BondedParamsAndLists bonded = parameters.getBondedParamsAndLists(psf);

  CHECK(bonded.paramsSize == std::vector<int>{3, 1, 2, 1, 1, 0});
  CHECK(bonded.listsSize == std::vector<int>{3, 1, 2, 1, 1, 0});

  const std::vector<std::vector<float>> expectedParams = {
      {1.1f, 100.0f},
      {1.2f, 200.0f},
      {1.3f, 300.0f},
      {2.5f, 10.0f},
      {static_cast<float>(apo_test::Deg2Rad(120.0)), 50.0f},
      {static_cast<float>(apo_test::Deg2Rad(110.0)), 60.0f},
      {3.0f, 0.5f, 0.0f, -1.0f},
      {180.0f, 1.5f, 0.0f, 1.0f},
  };

  apo_test::CheckVectorsClose2D<float>(bonded.paramsVal, expectedParams,
                                             FLOAT_TOLERANCE);

  const std::vector<std::vector<int>> expectedLists = {
      {0, 1, 0, 13},
      {1, 2, 1, 13},
      {2, 3, 2, 13},
      {0, 2, 0, 13},
      {0, 1, 2, 0, 13, 13},
      {1, 2, 3, 1, 13, 13},
      {0, 1, 2, 3, 0, 13, 13, 13},
      {0, 1, 2, 3, 0, 13, 13, 13},
  };

  CHECK(bonded.listVal == expectedLists);

  apo_test::RemoveIfExists(parameterFile);
  apo_test::RemoveIfExists(psfFile);
}

TEST_CASE("CharmmParametersBuildsVdwParamsAndTypes") {
  const std::string parameterFile = "tmp_charmm_parameters_full.prm";
  const std::string psfFile = "tmp_charmm_parameters_test.psf";
  apo_test::WriteTextFile(parameterFile, TEST_PARAMETER_TEXT);
  apo_test::WriteTextFile(psfFile, TEST_PSF_TEXT);

  CharmmParameters parameters(parameterFile);
  std::shared_ptr<CharmmPSF> psf = std::make_shared<CharmmPSF>(psfFile);

  const VdwParamsAndTypes vdw = parameters.getVdwParamsAndTypes(psf);

  CHECK(vdw.vdwTypes == std::vector<int>{0, 1, 2, 3});
  CHECK(vdw.vdw14Types == std::vector<int>{0, 1, 2, 3});

  // Four unique atom types produce a lower-triangular 4x4 packed matrix.
  REQUIRE(vdw.vdwParams.size() == 20);
  REQUIRE(vdw.vdw14Params.size() == 20);

  const double aaEpsilon = std::sqrt((-0.1) * (-0.1));
  const double aaRmin = 1.0 + 1.0;
  const float expectedAAC6 =
      static_cast<float>(2.0 * aaEpsilon * std::pow(aaRmin, 6));
  const float expectedAAC12 =
      static_cast<float>(aaEpsilon * std::pow(aaRmin, 12));

  CHECK(vdw.vdwParams[0] == Approx(expectedAAC6).margin(FLOAT_TOLERANCE));
  CHECK(vdw.vdwParams[1] == Approx(expectedAAC12).margin(FLOAT_TOLERANCE));

  // Pair A-D is packed at triangular pair index 6: (A,A), (B,A), (B,B),
  // (C,A), (C,B), (C,C), (D,A). It should use the NBFIX record.
  const std::size_t adPairIndex = 6;
  const double adEpsilon = 0.05;
  const double adRmin = 2.5;
  const float expectedADC6 =
      static_cast<float>(2.0 * adEpsilon * std::pow(adRmin, 6));
  const float expectedADC12 =
      static_cast<float>(adEpsilon * std::pow(adRmin, 12));

  CHECK(vdw.vdwParams[2 * adPairIndex + 0] ==
        Approx(expectedADC6).margin(FLOAT_TOLERANCE));
  CHECK(vdw.vdwParams[2 * adPairIndex + 1] ==
        Approx(expectedADC12).margin(FLOAT_TOLERANCE));

  // No A-D 1-4-specific NBFIX values were supplied, so the current
  // implementation uses the same NBFIX values for vdw14Params.
  CHECK(vdw.vdw14Params[2 * adPairIndex + 0] ==
        Approx(expectedADC6).margin(FLOAT_TOLERANCE));
  CHECK(vdw.vdw14Params[2 * adPairIndex + 1] ==
        Approx(expectedADC12).margin(FLOAT_TOLERANCE));

  apo_test::RemoveIfExists(parameterFile);
  apo_test::RemoveIfExists(psfFile);
}

TEST_CASE("CharmmParametersRejectsBadInputs") {
  SECTION("MissingFile") {
    const std::string parameterFile = "tmp_charmm_parameters_missing.prm";
    apo_test::RemoveIfExists(parameterFile);

    CHECK_THROWS_AS(CharmmParameters(parameterFile), std::invalid_argument);
  }

  SECTION("MalformedNonbondedLine") {
    const std::string parameterFile = "tmp_charmm_parameters_malformed.prm";
    apo_test::WriteTextFile(parameterFile, MALFORMED_PARAMETER_TEXT);

    CHECK_THROWS_AS(CharmmParameters(parameterFile), std::invalid_argument);

    apo_test::RemoveIfExists(parameterFile);
  }

  SECTION("MissingBondParameter") {
    const std::string parameterFile = "tmp_charmm_parameters_missing_bond.prm";
    const std::string psfFile = "tmp_charmm_parameters_test.psf";
    apo_test::WriteTextFile(parameterFile, MISSING_BOND_PARAMETER_TEXT);
    apo_test::WriteTextFile(psfFile, TEST_PSF_TEXT);

    CharmmParameters parameters(parameterFile);
    std::shared_ptr<CharmmPSF> psf = std::make_shared<CharmmPSF>(psfFile);

    CHECK_THROWS_AS(parameters.getBondedParamsAndLists(psf),
                    std::invalid_argument);

    apo_test::RemoveIfExists(parameterFile);
    apo_test::RemoveIfExists(psfFile);
  }
}
