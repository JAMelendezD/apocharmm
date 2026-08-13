// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apo_test_helpers.h"
#include "apocharmm_c/CharmmParameters.h"
#include "apocharmm_c/Error.h"
#include "apocharmm_c/Status.h"
#include "apocharmm_c/detail/CharmmParametersHandle.h"
#include "catch.hpp"

#include <cstddef>
#include <memory>
#include <string>

namespace {

const std::string VALID_PARAMETER_TEXT = "NONBONDED\n"
                                         "A 0.0 -0.1 1.0\n"
                                         "END\n";

const std::string SUPPLEMENTAL_PARAMETER_TEXT = "BONDS\n"
                                                "A B 100.0 1.0\n"
                                                "END\n";

const std::string MALFORMED_PARAMETER_TEXT = "NONBONDED\n"
                                             "A 0.0 -0.1 1.0 2.0\n"
                                             "END\n";

struct CharmmParametersDeleter {
  void operator()(apo_charmm_parameters *parameters) const noexcept {
    apo_charmm_parameters_destroy(parameters);
    return;
  }
};

using CharmmParametersHandle =
    std::unique_ptr<apo_charmm_parameters, CharmmParametersDeleter>;

CharmmParametersHandle MakeParameters(const std::string &parameterFile) {
  apo_charmm_parameters *parameters = nullptr;

  REQUIRE(apo_charmm_parameters_create(&parameters, parameterFile.c_str()) ==
          APO_STATUS_OK);
  REQUIRE(parameters != nullptr);

  return CharmmParametersHandle(parameters);
}

} // namespace

TEST_CASE("CapiCharmmParametersCreatesFromSingleAndMultipleFiles") {
  const std::string parameterFile = "tmp_capi_charmm_parameters_valid.prm";
  const std::string supplementalFile =
      "tmp_capi_charmm_parameters_supplemental.prm";

  apo_test::WriteTextFile(parameterFile, VALID_PARAMETER_TEXT);
  apo_test::WriteTextFile(supplementalFile, SUPPLEMENTAL_PARAMETER_TEXT);

  SECTION("SingleFile") {
    apo_charmm_parameters *rawParameters = nullptr;

    CHECK(apo_charmm_parameters_create(&rawParameters, parameterFile.c_str()) ==
          APO_STATUS_OK);
    REQUIRE(rawParameters != nullptr);

    CharmmParametersHandle parameters(rawParameters);
    CHECK(std::string(apo_last_error()).empty() == true);
  }

  SECTION("MultipleFiles") {
    const char *const paths[] = {parameterFile.c_str(),
                                 supplementalFile.c_str()};
    apo_charmm_parameters *rawParameters = nullptr;

    CHECK(apo_charmm_parameters_create_from_files(
              &rawParameters, paths, static_cast<std::size_t>(2)) ==
          APO_STATUS_OK);
    REQUIRE(rawParameters != nullptr);

    CharmmParametersHandle parameters(rawParameters);
    CHECK(std::string(apo_last_error()).empty() == true);
  }

  apo_test::RemoveIfExists(parameterFile);
  apo_test::RemoveIfExists(supplementalFile);
}

TEST_CASE("CapiCharmmParametersValidatesPointersCountsAndStrings") {
  const std::string parameterFile = "tmp_capi_charmm_parameters_validation.prm";
  apo_test::WriteTextFile(parameterFile, VALID_PARAMETER_TEXT);

  SECTION("SingleFile") {
    apo_test::CheckStatusAndDiagnostic(
        apo_charmm_parameters_create(nullptr, parameterFile.c_str()),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_charmm_parameters_create: out pointer is NULL");

    apo_charmm_parameters staleParameters;
    apo_charmm_parameters *parameters = &staleParameters;

    apo_test::CheckStatusAndDiagnostic(
        apo_charmm_parameters_create(&parameters, nullptr),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_charmm_parameters_create: path is NULL or empty");
    CHECK(parameters == nullptr);

    parameters = &staleParameters;

    apo_test::CheckStatusAndDiagnostic(
        apo_charmm_parameters_create(&parameters, ""),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_charmm_parameters_create: path is NULL or empty");
    CHECK(parameters == nullptr);
  }

  SECTION("MultipleFiles") {
    const char *const validPaths[] = {parameterFile.c_str()};

    apo_test::CheckStatusAndDiagnostic(
        apo_charmm_parameters_create_from_files(nullptr, validPaths, 1),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_charmm_parameters_create_from_files: out pointer is NULL");

    apo_charmm_parameters staleParameters;
    apo_charmm_parameters *parameters = &staleParameters;

    apo_test::CheckStatusAndDiagnostic(
        apo_charmm_parameters_create_from_files(&parameters, nullptr, 0),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_charmm_parameters_create_from_files: paths must contain at least "
        "one parameter file");
    CHECK(parameters == nullptr);

    parameters = &staleParameters;

    apo_test::CheckStatusAndDiagnostic(
        apo_charmm_parameters_create_from_files(&parameters, nullptr, 1),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_charmm_parameters_create_from_files: paths is NULL");
    CHECK(parameters == nullptr);

    const char *const nullPaths[] = {nullptr};
    parameters = &staleParameters;

    apo_test::CheckStatusAndDiagnostic(
        apo_charmm_parameters_create_from_files(&parameters, nullPaths, 1),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_charmm_parameters_create_from_files: path is NULL or empty");
    CHECK(parameters == nullptr);

    const char *const emptyPaths[] = {""};
    parameters = &staleParameters;

    apo_test::CheckStatusAndDiagnostic(
        apo_charmm_parameters_create_from_files(&parameters, emptyPaths, 1),
        APO_STATUS_INVALID_ARGUMENT,
        "apo_charmm_parameters_create_from_files: path is NULL or empty");
    CHECK(parameters == nullptr);
  }

  apo_test::RemoveIfExists(parameterFile);
}

TEST_CASE("CapiCharmmParametersMapsNativeErrors") {
  SECTION("MissingSingleFile") {
    const std::string parameterFile = "tmp_capi_charmm_parameters_missing.prm";
    apo_test::RemoveIfExists(parameterFile);

    apo_charmm_parameters staleParameters;
    apo_charmm_parameters *parameters = &staleParameters;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apo_charmm_parameters_create(
                       &parameters, parameterFile.c_str())));

    apo_test::CheckNativeError(
        status, APO_STATUS_RUNTIME_ERROR, "Runtime",
        "apo_charmm_parameters_create",
        "Failed to open CHARMM parameter file \"" + parameterFile + "\"",
        "src/CharmmParameters.cu", "readCharmmParameterFile");
    CHECK(parameters == nullptr);
  }

  SECTION("MalformedSingleFile") {
    const std::string parameterFile =
        "tmp_capi_charmm_parameters_malformed.prm";
    apo_test::WriteTextFile(parameterFile, MALFORMED_PARAMETER_TEXT);

    apo_charmm_parameters staleParameters;
    apo_charmm_parameters *parameters = &staleParameters;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apo_charmm_parameters_create(
                       &parameters, parameterFile.c_str())));

    apo_test::CheckNativeError(
        status, APO_STATUS_RUNTIME_ERROR, "Runtime",
        "apo_charmm_parameters_create",
        "Invalid NONBONDED parameter record in file \"" + parameterFile +
            "\" at line 2: A 0.0 -0.1 1.0 2.0",
        "src/CharmmParameters.cu", "parseNonbondedRecord");
    CHECK(parameters == nullptr);

    apo_test::RemoveIfExists(parameterFile);
  }

  SECTION("MissingFileInSequence") {
    const std::string parameterFile =
        "tmp_capi_charmm_parameters_valid_sequence.prm";
    const std::string missingFile =
        "tmp_capi_charmm_parameters_missing_sequence.prm";

    apo_test::WriteTextFile(parameterFile, VALID_PARAMETER_TEXT);
    apo_test::RemoveIfExists(missingFile);

    const char *const paths[] = {parameterFile.c_str(), missingFile.c_str()};
    apo_charmm_parameters staleParameters;
    apo_charmm_parameters *parameters = &staleParameters;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apo_charmm_parameters_create_from_files(
                       &parameters, paths, static_cast<std::size_t>(2))));

    apo_test::CheckNativeError(
        status, APO_STATUS_RUNTIME_ERROR, "Runtime",
        "apo_charmm_parameters_create_from_files",
        "Failed to open CHARMM parameter file \"" + missingFile + "\"",
        "src/CharmmParameters.cu", "readCharmmParameterFile");
    CHECK(parameters == nullptr);

    apo_test::RemoveIfExists(parameterFile);
  }
}

TEST_CASE("CapiCharmmParametersSuccessClearsStaleDiagnostic") {
  const std::string parameterFile =
      "tmp_capi_charmm_parameters_clear_error.prm";
  apo_test::WriteTextFile(parameterFile, VALID_PARAMETER_TEXT);

  REQUIRE(apo_charmm_parameters_create(nullptr, parameterFile.c_str()) ==
          APO_STATUS_INVALID_ARGUMENT);
  REQUIRE(std::string(apo_last_error()).empty() == false);

  CharmmParametersHandle parameters = MakeParameters(parameterFile);
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_test::RemoveIfExists(parameterFile);
}

TEST_CASE("CapiCharmmParametersDestroyIsNoexcept") {
  const std::string parameterFile = "tmp_capi_charmm_parameters_destroy.prm";
  apo_test::WriteTextFile(parameterFile, VALID_PARAMETER_TEXT);

  CharmmParametersHandle parameters = MakeParameters(parameterFile);
  apo_charmm_parameters *const rawParameters = parameters.release();

  CHECK_NOTHROW(apo_charmm_parameters_destroy(rawParameters));
  CHECK(std::string(apo_last_error()).empty() == true);

  CHECK_NOTHROW(apo_charmm_parameters_destroy(nullptr));
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_test::RemoveIfExists(parameterFile);
}
