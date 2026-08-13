// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apo_test_helpers.h"
#include "apocharmm_c/CharmmCrd.h"
#include "apocharmm_c/Error.h"
#include "apocharmm_c/Status.h"
#include "apocharmm_c/detail/CharmmCrdHandle.h"
#include "catch.hpp"

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

namespace {

struct CharmmCrdDeleter {
  void operator()(apo_charmm_crd *crd) const noexcept {
    apo_charmm_crd_destroy(crd);
    return;
  }
};

using CharmmCrdHandle = std::unique_ptr<apo_charmm_crd, CharmmCrdDeleter>;

class NegativeCountCharmmCrd : public CharmmCrd {
public:
  explicit NegativeCountCharmmCrd(const std::string &fileName)
      : CharmmCrd(fileName) {
    m_NumAtoms = -1;
  }
};

void WriteValidCharmmCrd(const std::string &fileName) {
  std::ofstream file(fileName);
  REQUIRE(file.good());

  file << "* generated C API CharmmCrd unit test\n";
  file << "*\n";
  file << std::setw(10) << 2 << "\n";

  file << std::setw(5) << 1 << std::setw(5) << 1 << " " << std::left
       << std::setw(4) << "RS1" << " " << std::setw(4) << "A1" << std::right
       << std::fixed << std::setprecision(5) << std::setw(10) << 1.25
       << std::setw(10) << -2.5 << std::setw(10) << 3.75 << "\n";

  file << std::setw(5) << 2 << std::setw(5) << 1 << " " << std::left
       << std::setw(4) << "RS1" << " " << std::setw(4) << "A2" << std::right
       << std::fixed << std::setprecision(5) << std::setw(10) << 0.0
       << std::setw(10) << 4.5 << std::setw(10) << -6.25 << "\n";

  file.close();
  REQUIRE(file.good());

  return;
}

CharmmCrdHandle MakeCrd(const std::string &crdFile) {
  apo_charmm_crd *crd = nullptr;

  REQUIRE(apo_charmm_crd_create(&crd, crdFile.c_str()) == APO_STATUS_OK);
  REQUIRE(crd != nullptr);

  return CharmmCrdHandle(crd);
}

} // namespace

TEST_CASE("CapiCharmmCrdCreatesAndReadsValues") {
  const std::string crdFile = "tmp_capi_charmm_crd_valid.cor";
  WriteValidCharmmCrd(crdFile);

  CharmmCrdHandle crd = MakeCrd(crdFile);

  std::size_t numAtoms = 0;
  CHECK(apo_charmm_crd_get_num_atoms(&numAtoms, crd.get()) == APO_STATUS_OK);
  CHECK(numAtoms == 2);

  std::vector<double> coordinates(3 * numAtoms);
  CHECK(apo_charmm_crd_get_coordinates(coordinates.data(), coordinates.size(),
                                       crd.get()) == APO_STATUS_OK);

  const std::vector<double> expected = {1.25, -2.5, 3.75, 0.0, 4.5, -6.25};
  apo_test::CheckVectorsClose<double>(coordinates, expected, 1.0e-12);
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_test::RemoveIfExists(crdFile);
}

TEST_CASE("CapiCharmmCrdValidatesCreationArguments") {
  const std::string crdFile = "tmp_capi_charmm_crd_create_validation.cor";
  WriteValidCharmmCrd(crdFile);

  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_crd_create(nullptr, crdFile.c_str()),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_crd_create: out pointer is NULL");

  apo_charmm_crd staleCrd;
  apo_charmm_crd *crd = &staleCrd;

  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_crd_create(&crd, nullptr), APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_crd_create: path is NULL or empty");
  CHECK(crd == nullptr);

  crd = &staleCrd;

  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_crd_create(&crd, ""), APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_crd_create: path is NULL or empty");
  CHECK(crd == nullptr);

  apo_test::RemoveIfExists(crdFile);
}

TEST_CASE("CapiCharmmCrdValidatesCountGetter") {
  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_crd_get_num_atoms(nullptr, nullptr),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_crd_get_num_atoms: num_atoms is NULL");

  std::size_t numAtoms = 99;
  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_crd_get_num_atoms(&numAtoms, nullptr),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_crd_get_num_atoms: CharmmCrd is NULL");
  CHECK(numAtoms == 0);

  apo_charmm_crd nullObjectCrd;
  numAtoms = 99;
  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_crd_get_num_atoms(&numAtoms, &nullObjectCrd),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_crd_get_num_atoms: CharmmCrd object is NULL");
  CHECK(numAtoms == 0);

  const std::string crdFile = "tmp_capi_charmm_crd_negative_count.cor";
  WriteValidCharmmCrd(crdFile);

  apo_charmm_crd negativeCountCrd;
  negativeCountCrd.object = std::make_shared<NegativeCountCharmmCrd>(crdFile);
  numAtoms = 99;
  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_crd_get_num_atoms(&numAtoms, &negativeCountCrd),
      APO_STATUS_RUNTIME_ERROR,
      "apo_charmm_crd_get_num_atoms: CharmmCrd returned a negative atom count");
  CHECK(numAtoms == 0);

  apo_test::RemoveIfExists(crdFile);
}

TEST_CASE("CapiCharmmCrdValidatesCoordinateBuffer") {
  const std::string crdFile = "tmp_capi_charmm_crd_buffer_validation.cor";
  WriteValidCharmmCrd(crdFile);
  CharmmCrdHandle crd = MakeCrd(crdFile);

  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_crd_get_coordinates(nullptr, 0, nullptr),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_crd_get_coordinates: CharmmCrd is NULL");

  apo_charmm_crd nullObjectCrd;
  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_crd_get_coordinates(nullptr, 0, &nullObjectCrd),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_crd_get_coordinates: CharmmCrd object is NULL");

  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_crd_get_coordinates(nullptr, 6, crd.get()),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_crd_get_coordinates: Coordinate output buffer is NULL");

  std::vector<double> coordinates(5);
  apo_test::CheckStatusAndDiagnostic(
      apo_charmm_crd_get_coordinates(coordinates.data(), coordinates.size(),
                                     crd.get()),
      APO_STATUS_INVALID_ARGUMENT,
      "apo_charmm_crd_get_coordinates: Coordinate output buffer is too small");

  apo_test::RemoveIfExists(crdFile);
}

TEST_CASE("CapiCharmmCrdMapsNativeErrors") {
  SECTION("MissingFile") {
    const std::string crdFile = "tmp_capi_charmm_crd_missing.cor";
    apo_test::RemoveIfExists(crdFile);

    apo_charmm_crd staleCrd;
    apo_charmm_crd *crd = &staleCrd;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apo_charmm_crd_create(&crd, crdFile.c_str())));

    apo_test::CheckNativeError(status, APO_STATUS_RUNTIME_ERROR, "Runtime",
                               "apo_charmm_crd_create",
                               "Failed to open file \"" + crdFile + "\"",
                               "src/str_utils.cpp", "read_file_into_string");
    CHECK(crd == nullptr);
  }

  SECTION("InvalidAtomCount") {
    const std::string crdFile = "tmp_capi_charmm_crd_invalid_count.cor";
    apo_test::WriteTextFile(
        crdFile,
        "* generated invalid-count C API CharmmCrd unit test\n*\nBAD\n");

    apo_charmm_crd staleCrd;
    apo_charmm_crd *crd = &staleCrd;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apo_charmm_crd_create(&crd, crdFile.c_str())));

    apo_test::CheckNativeError(
        status, APO_STATUS_RUNTIME_ERROR, "Runtime", "apo_charmm_crd_create",
        "Invalid atom count value \"BAD\" in CHARMM coordinate file \"" +
            crdFile + "\" at line 3",
        "src/CharmmCrd.cu", "readCharmmCrdFile");
    CHECK(crd == nullptr);

    apo_test::RemoveIfExists(crdFile);
  }

  SECTION("MalformedFile") {
    const std::string crdFile = "tmp_capi_charmm_crd_malformed.cor";
    apo_test::WriteTextFile(crdFile, "* generated malformed C API CharmmCrd "
                                     "unit test\n*\n         1\nshort\n");

    apo_charmm_crd staleCrd;
    apo_charmm_crd *crd = &staleCrd;
    apo_status status = APO_STATUS_OK;

    CHECK_NOTHROW((status = apo_charmm_crd_create(&crd, crdFile.c_str())));

    apo_test::CheckNativeError(
        status, APO_STATUS_RUNTIME_ERROR, "Runtime", "apo_charmm_crd_create",
        "Coordinate record 1 is truncated in CHARMM coordinate file \"" +
            crdFile + "\" at line 4",
        "src/CharmmCrd.cu", "ParseCharmmCrdCoordinate");
    CHECK(crd == nullptr);

    apo_test::RemoveIfExists(crdFile);
  }
}

TEST_CASE("CapiCharmmCrdSuccessClearsStaleDiagnostic") {
  const std::string crdFile = "tmp_capi_charmm_crd_clear_error.cor";
  WriteValidCharmmCrd(crdFile);

  REQUIRE(apo_charmm_crd_create(nullptr, crdFile.c_str()) ==
          APO_STATUS_INVALID_ARGUMENT);
  REQUIRE(std::string(apo_last_error()).empty() == false);

  CharmmCrdHandle crd = MakeCrd(crdFile);
  CHECK(std::string(apo_last_error()).empty() == true);

  std::size_t numAtoms = 0;
  CHECK(apo_charmm_crd_get_num_atoms(&numAtoms, crd.get()) == APO_STATUS_OK);
  CHECK(numAtoms == 2);
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_test::RemoveIfExists(crdFile);
}

TEST_CASE("CapiCharmmCrdDestroyIsNoexcept") {
  const std::string crdFile = "tmp_capi_charmm_crd_destroy.cor";
  WriteValidCharmmCrd(crdFile);

  CharmmCrdHandle crd = MakeCrd(crdFile);
  apo_charmm_crd *const rawCrd = crd.release();

  CHECK_NOTHROW(apo_charmm_crd_destroy(rawCrd));
  CHECK(std::string(apo_last_error()).empty() == true);

  CHECK_NOTHROW(apo_charmm_crd_destroy(nullptr));
  CHECK(std::string(apo_last_error()).empty() == true);

  apo_test::RemoveIfExists(crdFile);
}
