// BEGINLICENSE
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: James E. Gonzales II
//
// ENDLICENSE

#include "apocharmm_c/AtomSelection.h"
#include "apocharmm_c/detail/AtomSelectionHandle.h"
#include "apocharmm_c/detail/ErrorInternal.h"
#include "apocharmm_c/detail/Validation.h"

#include <string>
#include <vector>

extern "C" void apo_atom_selection_destroy(apo_atom_selection *selection) {
  const char *function_name = "apo_atom_selection_destroy";
  apocharmm_c::guard_destroy(
      [selection](void) -> void {
        delete selection;
        return;
      },
      function_name);
  return;
}

extern "C" apo_status
apo_atom_selection_get_num_atoms(size_t *num_atoms,
                                 const apo_atom_selection *selection) {
  const char *function_name = "apo_atom_selection_get_num_atoms";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<size_t>(
            num_atoms, function_name, "num_atoms"));

        *num_atoms = 0;

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_atom_selection>(
                selection, function_name, "AtomSelection"));

        const int n = selection->object->getNumAtoms();

        if (n < 0) {
          return apocharmm_c::set_last_error(
              APO_STATUS_RUNTIME_ERROR, function_name,
              "AtomSelection returned a negative atom count");
        }

        *num_atoms = static_cast<size_t>(n);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_atom_selection_get_num_selected(size_t *num_selected,
                                    const apo_atom_selection *selection) {
  const char *function_name = "apo_atom_selection_get_num_selected";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<size_t>(
            num_selected, function_name, "num_selected"));

        *num_selected = 0;

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_atom_selection>(
                selection, function_name, "AtomSelection"));

        const int n = selection->object->getNumSelected();

        if (n < 0) {
          return apocharmm_c::set_last_error(
              APO_STATUS_RUNTIME_ERROR,
              "apo_atom_selection_get_num_selected: AtomSelection returned a "
              "negative selected atom count");
        }

        *num_selected = static_cast<size_t>(n);

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_atom_selection_get_atom_indices(int *atom_indices, const size_t len,
                                    const apo_atom_selection *selection) {
  const char *function_name = "apo_atom_selection_get_atom_indices";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_atom_selection>(
                selection, function_name, "AtomSelection"));

        const std::vector<int> indices = selection->object->getAtomIndices();
        const size_t req_len = indices.size();

        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_output_buffer<int>(
            atom_indices, len, req_len, function_name, "Atom index buffer"));

        for (size_t i = 0; i < req_len; i++)
          atom_indices[i] = indices[i];

        return APO_STATUS_OK;
      },
      function_name);
}

extern "C" apo_status
apo_atom_selection_contains(bool *is_selected,
                            const apo_atom_selection *selection,
                            const int atom_index) {
  const char *function_name = "apo_atom_selection_contains";

  return apocharmm_c::guard(
      [&](void) -> apo_status {
        APOCHARMM_C_RETURN_IF_ERROR(apocharmm_c::require_pointer<bool>(
            is_selected, function_name, "is_selected"));

        *is_selected = false;

        APOCHARMM_C_RETURN_IF_ERROR(
            apocharmm_c::require_handle_object<apo_atom_selection>(
                selection, function_name, "AtomSelection"));

        const int num_atoms = selection->object->getNumAtoms();

        if ((atom_index < 0) || (atom_index >= num_atoms)) {
          return apocharmm_c::invalid_argument(function_name,
                                               "atom_index is out of range");
        }

        *is_selected = selection->object->contains(atom_index);

        return APO_STATUS_OK;
      },
      function_name);
}
