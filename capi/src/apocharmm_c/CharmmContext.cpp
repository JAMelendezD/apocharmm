// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apocharmm_c/CharmmContext.h"
#include "apocharmm_c/detail/CharmmContextHandle.h"
#include "apocharmm_c/detail/CharmmCrdHandle.h"
#include "apocharmm_c/detail/EnumConversion.h"
#include "apocharmm_c/detail/ErrorInternal.h"
#include "apocharmm_c/detail/ForceManagerHandle.h"
#include "apocharmm_c/detail/Validation.h"

#include <memory>

extern "C" apo_status
apo_charmm_context_create(apo_charmm_context **out,
                          const apo_force_manager *force_manager) {
  const char *function_name = "apo_charmm_context_create";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::prepare_output_pointer<apo_charmm_context>(
                out, function_name, "out"));

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_force_manager>(
                force_manager, function_name, "ForceManager"));

        std::unique_ptr<apo_charmm_context> handle(new apo_charmm_context());
        handle->force_manager = force_manager->object;
        handle->object = std::make_shared<CharmmContext>(handle->force_manager);

        *out = handle.release();

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" void apo_charmm_context_destroy(apo_charmm_context *context) {
  delete context;
  return;
}

extern "C" apo_status
apo_charmm_context_set_coordinates(apo_charmm_context *context,
                                   const apo_charmm_crd *crd) {
  const char *function_name = "apo_charmm_context_set_coordinates";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_crd>(
                crd, function_name, "CharmmCrd"));

        context->object->setCoordinates(crd->object);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_set_periodic_boundary_condition(apo_charmm_context *context,
                                                   const apo_pbc pbc) {
  const char *function_name =
      "apo_charmm_context_set_periodic_boundary_condition";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        PBC cpp_pbc = PBC::NONE;
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::to_pbc(&cpp_pbc, pbc, function_name));

        context->object->setPeriodicBoundaryCondition(cpp_pbc);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_set_random_seed_for_velocities(apo_charmm_context *context,
                                                  const uint64_t seed) {
  const char *function_name =
      "apo_charmm_context_set_random_seed_for_velocities";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        context->object->setRandomSeedForVelocities(seed);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_charmm_context_use_holonomic_constraints(
    apo_charmm_context *context, const bool useHolonomicConstraints) {
  const char *function_name = "apo_charmm_context_use_holonomic_constraints";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        context->object->useHolonomicConstraints(useHolonomicConstraints);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_get_num_atoms(size_t *num_atoms,
                                 const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_num_atoms";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<size_t>(
            num_atoms, function_name, "num_atoms"));

        *num_atoms = 0;

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        const int n = context->object->getNumAtoms();

        if (n < 0) {
          return apocharmm_c::set_last_error(
              APO_STATUS_RUNTIME_ERROR,
              "apo_charmm_context_get_num_atoms: CharmmContext returned a "
              "negative atom count");
        }

        *num_atoms = static_cast<size_t>(n);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_get_coordinates_charges(double *xyzq, const size_t xyzq_len,
                                           const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_coordinates_charges";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        CudaContainer<double4> &coordinatesCharges =
            context->object->getCoordinatesCharges();
        const size_t num_atoms = coordinatesCharges.size();
        const size_t req_len = 4 * num_atoms;

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_output_buffer<double>(
            xyzq, xyzq_len, req_len, function_name,
            "Coordinate and charge output buffer"));

        coordinatesCharges.transferToHost();

        for (size_t i = 0; i < num_atoms; i++) {
          xyzq[i * 4 + 0] = coordinatesCharges[i].x;
          xyzq[i * 4 + 1] = coordinatesCharges[i].y;
          xyzq[i * 4 + 2] = coordinatesCharges[i].z;
          xyzq[i * 4 + 3] = coordinatesCharges[i].w;
        }

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_get_box_dimensions(double *x, double *y, double *z,
                                      const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_box_dimensions";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        std::vector<double> box_dims = context->object->getBoxDimensions();

        *x = box_dims[0];
        *y = box_dims[1];
        *z = box_dims[2];

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_charmm_context_get_periodic_boundary_condition(
    apo_pbc *pbc, const apo_charmm_context *context) {
  const char *function_name =
      "apo_charmm_context_get_periodic_boundary_condition";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_pointer<apo_pbc>(pbc, function_name, "pbc"));

        return apocharmm_c::from_pbc(
            pbc, context->object->getPeriodicBoundaryCondition(),
            function_name);
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_assign_velocities_at_temperature(apo_charmm_context *context,
                                                    const double temperature) {
  const char *function_name =
      "apo_charmm_context_assign_velocities_at_temperature";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        context->object->assignVelocitiesAtTemperature(
            static_cast<float>(temperature));

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_compute_temperature(double *temperature,
                                       apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_compute_temperature";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            temperature, function_name, "temperature"));

        const float t = context->object->computeTemperature();

        *temperature = static_cast<double>(t);

        return APO_STATUS_OK;
      },
      function_name);
}
