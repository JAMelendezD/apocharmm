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
#include "apocharmm_c/detail/CharmmParametersHandle.h"
#include "apocharmm_c/detail/CharmmPsfHandle.h"
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

extern "C" apo_status apo_charmm_context_create_from_psf_parameters(
    apo_charmm_context **out, const apo_charmm_psf *psf,
    const apo_charmm_parameters *parameters) {
  const char *function_name = "apo_charmm_context_create_from_psf_parameters";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::prepare_output_pointer<apo_charmm_context>(
                out, function_name, "out"));

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_psf>(
                psf, function_name, "CharmmPsf"));

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_parameters>(
                parameters, function_name, "CharmmParameters"));

        std::unique_ptr<apo_charmm_context> handle(new apo_charmm_context());
        handle->psf = psf->object;
        handle->parameters = parameters->object;
        handle->object =
            std::make_shared<CharmmContext>(handle->psf, handle->parameters);
        handle->force_manager = handle->object->getForceManager();

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
apo_charmm_context_set_box_dimensions(apo_charmm_context *context,
                                      const double *box_dimensions,
                                      const size_t box_dimensions_len) {
  const char *function_name = "apo_charmm_context_set_box_dimensions";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            box_dimensions, function_name, "box_dimensions"));

        if (box_dimensions_len != 3) {
          return apocharmm_c::invalid_argument(
              function_name, "box_dimensions must contain exactly 3 elements");
        }

        const std::vector<double> cpp_box_dimensions = {
            box_dimensions[0], box_dimensions[1], box_dimensions[2]};

        context->object->setBoxDimensions(cpp_box_dimensions);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_charmm_context_set_kappa(apo_charmm_context *context,
                                                   const double kappa) {
  const char *function_name = "apo_charmm_context_set_kappa";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        if (!std::isfinite(kappa) || kappa < 0.0) {
          return apocharmm_c::invalid_argument(
              function_name, "kappa must be finite and nonnegative");
        }

        context->object->setKappa(static_cast<float>(kappa));

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_charmm_context_set_cutoff(apo_charmm_context *context,
                                                    const double cutoff) {
  const char *function_name = "apo_charmm_context_set_cutoff";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        if (!std::isfinite(cutoff) || cutoff <= 0.0) {
          return apocharmm_c::invalid_argument(
              function_name, "cutoff must be finite and positive");
        }

        context->object->setCutoff(static_cast<float>(cutoff));

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_charmm_context_set_ctonnb(apo_charmm_context *context,
                                                    const double ctonnb) {
  const char *function_name = "apo_charmm_context_set_ctonnb";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        if (!std::isfinite(ctonnb) || ctonnb <= 0.0) {
          return apocharmm_c::invalid_argument(
              function_name, "ctonnb must be finite and positive");
        }

        context->object->setCtonnb(static_cast<float>(ctonnb));

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status apo_charmm_context_set_ctofnb(apo_charmm_context *context,
                                                    const double ctofnb) {
  const char *function_name = "apo_charmm_context_set_ctofnb";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        if (!std::isfinite(ctofnb) || ctofnb <= 0.0) {
          return apocharmm_c::invalid_argument(
              function_name, "ctofnb must be finite and positive");
        }

        context->object->setCtofnb(static_cast<float>(ctofnb));

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_set_fft_grid(apo_charmm_context *context, const int *grid,
                                const size_t grid_len) {
  const char *function_name = "apo_charmm_context_set_fft_grid";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_pointer<int>(grid, function_name, "grid"));

        if (grid_len != 3) {
          return apocharmm_c::invalid_argument(
              function_name, "grid must contain exactly 3 elements");
        }

        context->object->setFFTGrid(grid[0], grid[1], grid[2]);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_set_pme_spline_order(apo_charmm_context *context,
                                        const int order) {
  const char *function_name = "apo_charmm_context_set_pme_spline_order";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        if (order <= 0) {
          return apocharmm_c::invalid_argument(function_name,
                                               "order must be positive");
        }

        context->object->setPmeSplineOrder(order);

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
apo_charmm_context_set_vdw_type(apo_charmm_context *context,
                                const int vdw_type) {
  const char *function_name = "apo_charmm_context_set_vdw_type";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        if ((vdw_type < 1) || (vdw_type > 6)) {
          return apocharmm_c::invalid_argument(function_name,
                                               "vdw_type must be [1, 6]");
        }

        context->object->setVdwType(vdw_type);

        return APO_STATUS_OK;
      },
      function_name);
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
apo_charmm_context_get_velocity_mass(double *xyzm, const size_t xyzm_len,
                                     const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_velocity_mass";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        CudaContainer<double4> &velocityMass =
            context->object->getVelocityMass();
        const size_t num_atoms = velocityMass.size();
        const size_t req_len = 4 * num_atoms;

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_output_buffer<double>(
            xyzm, xyzm_len, req_len, function_name,
            "Velocity and mass output buffer"));

        velocityMass.transferToHost();

        for (size_t i = 0; i < num_atoms; i++) {
          xyzm[i * 4 + 0] = velocityMass[i].x;
          xyzm[i * 4 + 1] = velocityMass[i].y;
          xyzm[i * 4 + 2] = velocityMass[i].z;
          xyzm[i * 4 + 3] = velocityMass[i].w;
        }

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_get_box_dimensions(double *box_dimensions,
                                      const size_t box_dimensions_len,
                                      const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_box_dimensions";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        std::vector<double> box_dims = context->object->getBoxDimensions();

        if (box_dims.size() != 3) {
          return apocharmm_c::set_last_error(
              APO_STATUS_RUNTIME_ERROR,
              "apo_charmm_context_get_box_dimensions: CharmmContext did not "
              "return exactly 3 box dimensions");
        }

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_output_buffer<double>(
            box_dimensions, box_dimensions_len, 3, function_name,
            "Box dimension buffer"));

        for (size_t i = 0; i < 3; i++)
          box_dimensions[i] = box_dims[i];

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_get_kappa(double *kappa, const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_kappa";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            kappa, function_name, "kappa"));

        *kappa = static_cast<double>(context->object->getKappa());

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_get_cutoff(double *cutoff,
                              const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_cutoff";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            cutoff, function_name, "cutoff"));

        *cutoff = static_cast<double>(context->object->getCutoff());

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_get_ctonnb(double *ctonnb,
                              const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_ctonnb";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            ctonnb, function_name, "ctonnb"));

        *ctonnb = static_cast<double>(context->object->getCtonnb());

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_get_ctofnb(double *ctofnb,
                              const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_ctofnb";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<double>(
            ctofnb, function_name, "ctofnb"));

        *ctofnb = static_cast<double>(context->object->getCtofnb());

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_get_fft_grid(int *grid, const size_t grid_len,
                                const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_fft_grid";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        std::vector<int> fft_grid = context->object->getFFTGrid();

        if (fft_grid.size() != 3) {
          return apocharmm_c::set_last_error(
              APO_STATUS_RUNTIME_ERROR,
              "apo_charmm_context_get_fft_grid: CharmmContext did not return "
              "exactly 3 FFT grid dimensions");
        }

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_output_buffer<int>(
            grid, grid_len, 3, function_name, "FFT grid buffer"));

        for (size_t i = 0; i < 3; i++)
          grid[i] = fft_grid[i];

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_charmm_context_get_pme_spline_order(int *order,
                                        const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_pme_spline_order";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_pointer<int>(order, function_name, "order"));

        *order = context->object->getPmeSplineOrder();

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
apo_charmm_context_get_vdw_type(int *vdw_type,
                                const apo_charmm_context *context) {
  const char *function_name = "apo_charmm_context_get_vdw_type";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_charmm_context>(
                context, function_name, "CharmmContext"));

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<int>(
            vdw_type, function_name, "vdw_type"));

        *vdw_type = context->object->getVdwType();

        return APO_STATUS_OK;
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
