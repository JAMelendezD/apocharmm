// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apocharmm_c/HarmonicRestraintForce.h"
#include "apocharmm_c/detail/AtomSelectionHandle.h"
#include "apocharmm_c/detail/ErrorInternal.h"
#include "apocharmm_c/detail/ForceManagerHandle.h"
#include "apocharmm_c/detail/HarmonicRestraintForceHandle.h"
#include "apocharmm_c/detail/Validation.h"

#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

extern "C" apo_status
apo_harmonic_restraint_force_create(apo_harmonic_restraint_force **out,
                                    const int num_atoms) {
  const char *function_name = "apo_harmonic_restraint_force_create";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::prepare_output_pointer<apo_harmonic_restraint_force>(
                out, function_name, "out"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_positive_int(
            num_atoms, function_name, "num_atoms"));

        std::unique_ptr<apo_harmonic_restraint_force> handle(
            new apo_harmonic_restraint_force());

        handle->object =
            std::make_shared<HarmonicRestraintForce<long long int, float>>(
                num_atoms);

        *out = handle.release();

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" void
apo_harmonic_restraint_force_destroy(apo_harmonic_restraint_force *restraint) {
  delete restraint;
  return;
}

extern "C" apo_status apo_harmonic_restraint_force_set_selection(
    apo_harmonic_restraint_force *restraint,
    const apo_atom_selection *selection) {
  const char *function_name = "apo_harmonic_restraint_force_set_selection";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_harmonic_restraint_force>(
                restraint, function_name, "HarmonicRestraintForce"));

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_atom_selection>(
                selection, function_name, "AtomSelection"));

        restraint->object->setSelection(*selection->object);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_harmonic_restraint_force_set_force_constant(
    apo_harmonic_restraint_force *restraint, const double force_constant) {
  const char *function_name = "apo_harmonic_restraint_force_set_force_constant";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_harmonic_restraint_force>(
                restraint, function_name, "HarmonicRestraintForce"));

        if (!std::isfinite(force_constant)) {
          return apocharmm_c::invalid_argument(function_name,
                                               "force_constant must be finite");
        }

        if (force_constant < 0.0) {
          return apocharmm_c::invalid_argument(
              function_name, "force_constant must be non-negative");
        }

        restraint->object->setForceConstant(force_constant);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_harmonic_restraint_force_set_force_constants(
    apo_harmonic_restraint_force *restraint, const double *force_constants,
    const size_t force_constants_len) {
  const char *function_name =
      "apo_harmonic_restraint_force_set_force_constants";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_harmonic_restraint_force>(
                restraint, function_name, "HarmonicRestraintForce"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            force_constants, function_name, "force_constants"));

        std::vector<double> cpp_force_constants(force_constants_len);
        for (size_t i = 0; i < force_constants_len; i++)
          cpp_force_constants[i] = force_constants[i];

        restraint->object->setForceConstants(cpp_force_constants);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_harmonic_restraint_force_set_reference_coordinates(
    apo_harmonic_restraint_force *restraint,
    const double *reference_coordinates,
    const size_t reference_coordinates_len) {
  const char *function_name =
      "apo_harmonic_restraint_force_set_reference_coordinates";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_harmonic_restraint_force>(
                restraint, function_name, "HarmonicRestraintForce"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            reference_coordinates, function_name, "reference_coordinates"));

        if (reference_coordinates_len % 3 != 0) {
          return apocharmm_c::invalid_argument(
              function_name,
              "reference_coordinates length must be divisible by 3");
        }

        const size_t num_coordiantes = reference_coordinates_len / 3;

        std::vector<std::vector<double>> cpp_reference_coordinates(
            num_coordiantes, std::vector<double>(3));
        for (size_t i = 0; i < num_coordiantes; i++) {
          for (size_t j = 0; j < 3; j++)
            cpp_reference_coordinates[i][j] = reference_coordinates[i * 3 + j];
        }

        restraint->object->setReferenceCoordinates(cpp_reference_coordinates);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_harmonic_restraint_force_set_masses(apo_harmonic_restraint_force *restraint,
                                        const double *masses,
                                        const size_t masses_len) {
  const char *function_name = "apo_harmonic_restraint_force_set_masses";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_harmonic_restraint_force>(
                restraint, function_name, "HarmonicRestraintForce"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            masses, function_name, "masses"));

        std::vector<double> cpp_masses(masses_len);
        for (size_t i = 0; i < masses_len; i++)
          cpp_masses[i] = masses[i];

        restraint->object->setMasses(cpp_masses);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_harmonic_restraint_force_set_box_dimensions(
    apo_harmonic_restraint_force *restraint, const double *box_dimensions,
    const size_t box_dimensions_len) {
  const char *function_name = "apo_harmonic_restraint_force_set_box_dimensions";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_harmonic_restraint_force>(
                restraint, function_name, "HarmonicRestraintForce"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            box_dimensions, function_name, "box_dimensions"));

        if (box_dimensions_len != 3) {
          return apocharmm_c::invalid_argument(
              function_name, "box_dimensions must contain exactly 3 elements");
        }

        std::vector<double> cpp_box_dimensions(3);
        for (size_t i = 0; i < 3; i++) {
          if (!std::isfinite(box_dimensions[i])) {
            return apocharmm_c::invalid_argument(
                function_name, "box_dimensions values must be finite");
          }

          if (box_dimensions[i] <= 0.0) {
            return apocharmm_c::invalid_argument(
                function_name, "box_dimensions values must be positive");
          }

          cpp_box_dimensions[i] = box_dimensions[i];
        }

        restraint->object->setBoxDimensions(cpp_box_dimensions);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_force_manager_subscribe_harmonic_restraint_force(
    apo_force_manager *force_manager, apo_harmonic_restraint_force *restraint,
    const char *force_tag) {
  const char *function_name =
      "apo_force_manager_subscribe_harmonic_restraint_force";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_force_manager>(
                force_manager, function_name, "ForceManager"));

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_harmonic_restraint_force>(
                restraint, function_name, "HarmonicRestraintForce"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_c_string(
            force_tag, function_name, "force_tag"));

        force_manager->object->subscribe(
            restraint->object, std::string(force_tag),
            restraint->object->getStream(), restraint->object->getForce(),
            restraint->object->getEnergyVirial());

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_force_manager_unsubscribe_harmonic_restraint_force(
    apo_force_manager *force_manager, apo_harmonic_restraint_force *restraint) {
  const char *function_name =
      "apo_force_manager_unsubscribe_harmonic_restraint_force";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_force_manager>(
                force_manager, function_name, "ForceManager"));

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_harmonic_restraint_force>(
                restraint, function_name, "HarmonicRestraintForce"));

        force_manager->object->unsubscribe(restraint->object);

        return APO_STATUS_OK;
      },
      function_name);
}
