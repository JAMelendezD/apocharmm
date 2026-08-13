// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apo_test_helpers.h"
#include "apocharmm_c/CharmmPsf.h"
#include "apocharmm_c/Error.h"
#include "apocharmm_c/Status.h"
#include "apocharmm_c/detail/CharmmPsfHandle.h"
#include "catch.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace {

const std::string VALID_PSF_TEXT = R"PSF(PSF

       1 !NTITLE
 REMARKS generated C API CharmmPsf unit test
       2 !NATOM
       1 SEG1     1 RES1 C1   CT1   -0.500000  12.0110           0
       2 SEG1     1 RES1 H1   HA1    0.500000   1.0080           0
       1 !NBOND: bonds
       1       2
       0 !NTHETA: angles
       0 !NPHI: dihedrals
       0 !NIMPHI: impropers
       0 !NDON: donors
       0 !NACC: acceptors
       0 !NCRTERM: cross-terms
)PSF";

const std::string MALFORMED_PSF_TEXT = R"PSF(PSF

       1 !NTITLE
 REMARKS generated malformed C API CharmmPsf unit test
       1 !NATOM
       1 SEG1     1 RES1 C1   CT1    0.000000  12.0110           0
)PSF";

struct CharmmPsfDeleter {
  void operator()(apo_charmm_psf *psf) const noexcept {
    apo_charmm_psf_destroy(psf);
    return;
  }
};

using CharmmPsfHandle = std::unique_ptr<apo_charmm_psf, CharmmPsfDeleter>;
using CountGetter = apo_status (*)(std::size_t *, const apo_charmm_psf *);
using ScalarGetter = apo_status (*)(double *, const apo_charmm_psf *);

template <typename T>
using BufferGetter = apo_status (*)(T *, std::size_t, const apo_charmm_psf *);

CharmmPsfHandle MakePsf(const std::string &psfFile) {
  apo_charmm_psf *psf = nullptr;

  REQUIRE(apo_charmm_psf_create(&psf, psfFile.c_str()) == APO_STATUS_OK);
  REQUIRE(psf != nullptr);

  return CharmmPsfHandle(psf);
}

void CheckCountGetterValidation(const CountGetter getter,
                                const std::string &functionName,
                                const std::string &outputName,
                                const std::string &negativeCountMessage) {
  apo_test::CheckStatusAndDiagnostic(
      getter(nullptr, nullptr), APO_STATUS_INVALID_ARGUMENT,
      functionName + ": " + outputName + " is NULL");

  std::size_t value = 99;
  apo_test::CheckStatusAndDiagnostic(getter(&value, nullptr),
                                     APO_STATUS_INVALID_ARGUMENT,
                                     functionName + ": CharmmPsf is NULL");
  CHECK(value == 0);

  apo_charmm_psf nullObjectPsf;
  value = 99;
  apo_test::CheckStatusAndDiagnostic(
      getter(&value, &nullObjectPsf), APO_STATUS_INVALID_ARGUMENT,
      functionName + ": CharmmPsf object is NULL");
  CHECK(value == 0);

  apo_charmm_psf defaultPsf;
  defaultPsf.object = std::make_shared<CharmmPSF>();
  value = 99;
  apo_test::CheckStatusAndDiagnostic(getter(&value, &defaultPsf),
                                     APO_STATUS_RUNTIME_ERROR,
                                     negativeCountMessage);
  CHECK(value == 0);

  return;
}

template <typename T>
void CheckBufferGetterValidation(const BufferGetter<T> getter,
                                 const std::string &functionName,
                                 const std::string &bufferName,
                                 const std::size_t requiredLength,
                                 const apo_charmm_psf *psf) {
  REQUIRE(requiredLength > 0);

  apo_test::CheckStatusAndDiagnostic(getter(nullptr, 0, nullptr),
                                     APO_STATUS_INVALID_ARGUMENT,
                                     functionName + ": CharmmPsf is NULL");

  apo_charmm_psf nullObjectPsf;
  apo_test::CheckStatusAndDiagnostic(
      getter(nullptr, 0, &nullObjectPsf), APO_STATUS_INVALID_ARGUMENT,
      functionName + ": CharmmPsf object is NULL");

  apo_test::CheckStatusAndDiagnostic(
      getter(nullptr, requiredLength, psf), APO_STATUS_INVALID_ARGUMENT,
      functionName + ": " + bufferName + " is NULL");

  std::vector<T> buffer(requiredLength);
  apo_test::CheckStatusAndDiagnostic(
      getter(buffer.data(), requiredLength - 1, psf),
      APO_STATUS_INVALID_ARGUMENT,
      functionName + ": " + bufferName + " is too small");

  return;
}

void CheckScalarGetterValidation(const ScalarGetter getter,
                                 const std::string &functionName,
                                 const std::string &outputName) {
  apo_test::CheckStatusAndDiagnostic(
      getter(nullptr, nullptr), APO_STATUS_INVALID_ARGUMENT,
      functionName + ": " + outputName + " is NULL");

  double value = 99.0;
  apo_test::CheckStatusAndDiagnostic(getter(&value, nullptr),
                                     APO_STATUS_INVALID_ARGUMENT,
                                     functionName + ": CharmmPsf is NULL");
  CHECK(value == 0.0);

  apo_charmm_psf nullObjectPsf;
  value = 99.0;
  apo_test::CheckStatusAndDiagnostic(
      getter(&value, &nullObjectPsf), APO_STATUS_INVALID_ARGUMENT,
      functionName + ": CharmmPsf object is NULL");
  CHECK(value == 0.0);

  return;
}

} // namespace

TEST_CASE("CapiCharmmPsfCreatesAndReadsValues") {
  const std::string psfFile = "tmp_capi_charmm_psf_valid.psf";
  apo_test::WriteTextFile(psfFile, VALID_PSF_TEXT);

  CharmmPsfHandle psf = MakePsf(psfFile);

  std::size_t count = 0;
  CHECK(apo_charmm_psf_get_num_atoms(&count, psf.get()) == APO_STATUS_OK);
  CHECK(count == 2);
  CHECK(apo_charmm_psf_get_num_bonds(&count, psf.get()) == APO_STATUS_OK);
  CHECK(count == 1);
  CHECK(apo_charmm_psf_get_num_angles(&count, psf.get()) == APO_STATUS_OK);
  CHECK(count == 0);
  CHECK(apo_charmm_psf_get_num_dihedrals(&count, psf.get()) == APO_STATUS_OK);
  CHECK(count == 0);
  CHECK(apo_charmm_psf_get_num_impropers(&count, psf.get()) == APO_STATUS_OK);
  CHECK(count == 0);
  CHECK(apo_charmm_psf_get_num_cross_terms(&count, psf.get()) == APO_STATUS_OK);
  CHECK(count == 0);

  std::vector<char> stringBuffer(16);
  CHECK(apo_charmm_psf_get_segment_identifiers(stringBuffer.data(),
                                               stringBuffer.size(),
                                               psf.get()) == APO_STATUS_OK);
  CHECK(std::string(stringBuffer.data(), 8) == "SEG1    ");
  CHECK(std::string(stringBuffer.data() + 8, 8) == "SEG1    ");

  std::vector<int> residueIdentifiers(2);
  CHECK(apo_charmm_psf_get_residue_identifiers(residueIdentifiers.data(),
                                               residueIdentifiers.size(),
                                               psf.get()) == APO_STATUS_OK);
  const std::vector<int> expectedResidueIdentifiers = {1, 1};
  CHECK(residueIdentifiers == expectedResidueIdentifiers);

  CHECK(apo_charmm_psf_get_residue_names(stringBuffer.data(),
                                         stringBuffer.size(),
                                         psf.get()) == APO_STATUS_OK);
  CHECK(std::string(stringBuffer.data(), 8) == "RES1    ");
  CHECK(std::string(stringBuffer.data() + 8, 8) == "RES1    ");

  CHECK(apo_charmm_psf_get_atom_names(stringBuffer.data(), stringBuffer.size(),
                                      psf.get()) == APO_STATUS_OK);
  CHECK(std::string(stringBuffer.data(), 8) == "C1      ");
  CHECK(std::string(stringBuffer.data() + 8, 8) == "H1      ");

  CHECK(apo_charmm_psf_get_atom_types(stringBuffer.data(), stringBuffer.size(),
                                      psf.get()) == APO_STATUS_OK);
  CHECK(std::string(stringBuffer.data(), 8) == "CT1     ");
  CHECK(std::string(stringBuffer.data() + 8, 8) == "HA1     ");

  std::vector<double> values(2);
  CHECK(apo_charmm_psf_get_charges(values.data(), values.size(), psf.get()) ==
        APO_STATUS_OK);
  CHECK(values[0] == Approx(-0.5));
  CHECK(values[1] == Approx(0.5));

  CHECK(apo_charmm_psf_get_masses(values.data(), values.size(), psf.get()) ==
        APO_STATUS_OK);
  CHECK(values[0] == Approx(12.011));
  CHECK(values[1] == Approx(1.008));

  double scalar = 0.0;
  CHECK(apo_charmm_psf_get_net_charge(&scalar, psf.get()) == APO_STATUS_OK);
  CHECK(scalar == Approx(0.0));
  CHECK(apo_charmm_psf_get_total_mass(&scalar, psf.get()) == APO_STATUS_OK);
  CHECK(scalar == Approx(13.019));

  std::vector<char> fileNameBuffer(psfFile.size());
  CHECK(apo_charmm_psf_get_file_name(fileNameBuffer.data(),
                                     fileNameBuffer.size(),
                                     psf.get()) == APO_STATUS_OK);
  CHECK(std::string(fileNameBuffer.data(), fileNameBuffer.size()) == psfFile);
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_test::RemoveIfExists(psfFile);
}

TEST_CASE("CapiCharmmPsfValidatesCreationArguments") {
  const std::string psfFile = "tmp_capi_charmm_psf_create_validation.psf";
  apo_test::WriteTextFile(psfFile, VALID_PSF_TEXT);

  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_psf_create(nullptr, psfFile.c_str()),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_psf_create: out pointer is NULL");

  apo_charmm_psf stalePsf;
  apo_charmm_psf *psf = &stalePsf;
  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_psf_create(&psf, nullptr), APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_psf_create: PSF path is NULL or empty");
  CHECK(psf == nullptr);

  psf = &stalePsf;
  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_psf_create(&psf, ""), APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_psf_create: PSF path is NULL or empty");
  CHECK(psf == nullptr);

  apo_test::RemoveIfExists(psfFile);
}

TEST_CASE("CapiCharmmPsfValidatesCountGetters") {
  CheckCountGetterValidation(
      apo_charmm_psf_get_num_atoms, "apo_charmm_psf_get_num_atoms", "num_atoms",
      "apo_charmm_psf_get_num_atoms: CharmmPSF returned a negative "
      "atom count");
  CheckCountGetterValidation(
      apo_charmm_psf_get_num_bonds, "apo_charmm_psf_get_num_bonds", "num_bonds",
      "apo_charmm_psf_get_num_bonds: CharmmPSF returned a negative "
      "bond count");
  CheckCountGetterValidation(
      apo_charmm_psf_get_num_angles, "apo_charmm_psf_get_num_angles",
      "num_angles",
      "apo_charmm_psf_get_num_angles: CharmmPSF returned a negative "
      "angle count");
  CheckCountGetterValidation(
      apo_charmm_psf_get_num_dihedrals, "apo_charmm_psf_get_num_dihedrals",
      "num_dihedrals",
      "apo_charmm_psf_get_num_dihedrals: CharmmPSF returned a negative "
      "dihedral count");
  CheckCountGetterValidation(
      apo_charmm_psf_get_num_impropers, "apo_charmm_psf_get_num_impropers",
      "num_impropers",
      "apo_charmm_psf_get_num_impropers: CharmmPSF returned a negative "
      "improper count");
  CheckCountGetterValidation(
      apo_charmm_psf_get_num_cross_terms, "apo_charmm_psf_get_num_cross_terms",
      "num_cross_terms",
      "apo_charmm_psf_get_num_cross_terms: CharmmPSF returned a negative "
      "cross_term count");
}

TEST_CASE("CapiCharmmPsfValidatesBufferAndScalarGetters") {
  const std::string psfFile = "tmp_capi_charmm_psf_getter_validation.psf";
  apo_test::WriteTextFile(psfFile, VALID_PSF_TEXT);
  CharmmPsfHandle psf = MakePsf(psfFile);

  CheckBufferGetterValidation<char>(apo_charmm_psf_get_segment_identifiers,
                                    "apo_charmm_psf_get_segment_identifiers",
                                    "Segment identifier buffer", 16, psf.get());
  CheckBufferGetterValidation<int>(apo_charmm_psf_get_residue_identifiers,
                                   "apo_charmm_psf_get_residue_identifiers",
                                   "Residue identifier buffer", 2, psf.get());
  CheckBufferGetterValidation<char>(apo_charmm_psf_get_residue_names,
                                    "apo_charmm_psf_get_residue_names",
                                    "Residue name buffer", 16, psf.get());
  CheckBufferGetterValidation<char>(apo_charmm_psf_get_atom_names,
                                    "apo_charmm_psf_get_atom_names",
                                    "Atom name buffer", 16, psf.get());
  CheckBufferGetterValidation<char>(apo_charmm_psf_get_atom_types,
                                    "apo_charmm_psf_get_atom_types",
                                    "Atom type buffer", 16, psf.get());
  CheckBufferGetterValidation<double>(apo_charmm_psf_get_charges,
                                      "apo_charmm_psf_get_charges",
                                      "Charge buffer", 2, psf.get());
  CheckBufferGetterValidation<double>(apo_charmm_psf_get_masses,
                                      "apo_charmm_psf_get_masses",
                                      "Mass buffer", 2, psf.get());
  CheckBufferGetterValidation<char>(
      apo_charmm_psf_get_file_name, "apo_charmm_psf_get_file_name",
      "File name buffer", psfFile.size(), psf.get());

  CheckScalarGetterValidation(apo_charmm_psf_get_net_charge,
                              "apo_charmm_psf_get_net_charge", "net_charge");
  CheckScalarGetterValidation(apo_charmm_psf_get_total_mass,
                              "apo_charmm_psf_get_total_mass", "total_mass");

  apo_test::RemoveIfExists(psfFile);
}

TEST_CASE("CapiCharmmPsfMapsNativeErrors") {
  SECTION("MissingFile") {
    const std::string psfFile = "tmp_capi_charmm_psf_missing.psf";
    apo_test::RemoveIfExists(psfFile);

    apo_charmm_psf stalePsf;
    apo_charmm_psf *psf = &stalePsf;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apo_charmm_psf_create(&psf, psfFile.c_str())));

    apo_test::CheckNativeError(status, APO_STATUS_RUNTIME_ERROR, "Runtime",
                               "apo_charmm_psf_create",
                               "Failed to open file \"" + psfFile + "\"",
                               "src/str_utils.cpp", "read_file_into_string");
    CHECK(psf == nullptr);
  }

  SECTION("MalformedFile") {
    const std::string psfFile = "tmp_capi_charmm_psf_malformed.psf";
    apo_test::WriteTextFile(psfFile, MALFORMED_PSF_TEXT);

    apo_charmm_psf stalePsf;
    apo_charmm_psf *psf = &stalePsf;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apo_charmm_psf_create(&psf, psfFile.c_str())));

    apo_test::CheckNativeError(
        status, APO_STATUS_RUNTIME_ERROR, "Runtime", "apo_charmm_psf_create",
        "Could not find BOND section in PSF \"" + psfFile + "\"",
        "src/CharmmPSF.cu", "FindPsfSection");
    CHECK(psf == nullptr);

    apo_test::RemoveIfExists(psfFile);
  }

  SECTION("UninitializedAggregate") {
    apo_charmm_psf psf;
    psf.object = std::make_shared<CharmmPSF>();

    double value = 99.0;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apo_charmm_psf_get_net_charge(&value, &psf)));

    apo_test::CheckNativeError(status, APO_STATUS_NOT_INITIALIZED,
                               "NotInitialized",
                               "apo_charmm_psf_get_net_charge",
                               "CharmmPSF atom count is not initialized",
                               "src/CharmmPSF.cu", "getNetCharge");
    CHECK(value == 0.0);
  }

  SECTION("InconsistentChargeCount") {
    apo_charmm_psf psf;
    psf.object = std::make_shared<CharmmPSF>();
    psf.object->setNumAtoms(2);
    psf.object->getCharges().resize(1);

    double value = 99.0;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apo_charmm_psf_get_net_charge(&value, &psf)));

    apo_test::CheckNativeError(
        status, APO_STATUS_RUNTIME_ERROR, "Runtime",
        "apo_charmm_psf_get_net_charge",
        "CharmmPSF charge count does not match atom count; expected 2, "
        "observed 1",
        "src/CharmmPSF.cu", "getNetCharge");
    CHECK(value == 0.0);
  }

  SECTION("InconsistentMassCount") {
    apo_charmm_psf psf;
    psf.object = std::make_shared<CharmmPSF>();
    psf.object->setNumAtoms(2);
    psf.object->getMasses().resize(1);

    double value = 99.0;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apo_charmm_psf_get_total_mass(&value, &psf)));

    apo_test::CheckNativeError(
        status, APO_STATUS_RUNTIME_ERROR, "Runtime",
        "apo_charmm_psf_get_total_mass",
        "CharmmPSF mass count does not match atom count; expected 2, "
        "observed 1",
        "src/CharmmPSF.cu", "getTotalMass");
    CHECK(value == 0.0);
  }
}

TEST_CASE("CapiCharmmPsfSuccessClearsStaleDiagnostic") {
  const std::string psfFile = "tmp_capi_charmm_psf_clear_error.psf";
  apo_test::WriteTextFile(psfFile, VALID_PSF_TEXT);

  REQUIRE(apo_charmm_psf_create(nullptr, psfFile.c_str()) ==
          APO_STATUS_INVALID_ARGUMENT);
  REQUIRE(std::string(apo_last_error()).empty() == false);

  CharmmPsfHandle psf = MakePsf(psfFile);
  CHECK(std::string(apo_last_error()).empty() == true);

  std::size_t numAtoms = 0;
  CHECK(apo_charmm_psf_get_num_atoms(&numAtoms, psf.get()) == APO_STATUS_OK);
  CHECK(numAtoms == 2);
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_test::RemoveIfExists(psfFile);
}

TEST_CASE("CapiCharmmPsfDestroyIsNoexcept") {
  const std::string psfFile = "tmp_capi_charmm_psf_destroy.psf";
  apo_test::WriteTextFile(psfFile, VALID_PSF_TEXT);

  CharmmPsfHandle psf = MakePsf(psfFile);
  apo_charmm_psf *const rawPsf = psf.release();

  CHECK_NOTHROW(apo_charmm_psf_destroy(rawPsf));
  CHECK(std::string(apo_last_error()).empty() == true);

  CHECK_NOTHROW(apo_charmm_psf_destroy(nullptr));
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_test::RemoveIfExists(psfFile);
}
