# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

import ctypes

from ._base import _ApoObject
from ._lib import encode_path, lib
from ._types import FilePath
from .error import configure_status_function

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_charmm_psf_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_char_p],
        "CharmmPsf construction",
    )

    lib().apo_charmm_psf_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_charmm_psf_destroy.restype = None

    configure_status_function(
        lib().apo_charmm_psf_get_num_atoms,
        [ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p],
        "CharmmPsf.getNumAtoms()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_num_bonds,
        [ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p],
        "CharmmPsf.getNumBonds()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_num_angles,
        [ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p],
        "CharmmPsf.getNumAngles()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_num_dihedrals,
        [ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p],
        "CharmmPsf.getNumDihedrals()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_num_impropers,
        [ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p],
        "CharmmPsf.getNumImpropers()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_num_cross_terms,
        [ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p],
        "CharmmPsf.getNumCrossTerms()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_segment_identifiers,
        [ctypes.POINTER(ctypes.c_char), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmPsf.getSegmentIdentifiers()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_residue_identifiers,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmPsf.getResidueIdentifiers()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_residue_names,
        [ctypes.POINTER(ctypes.c_char), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmPsf.getResidueNames()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_atom_names,
        [ctypes.POINTER(ctypes.c_char), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmPsf.getAtomNames()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_atom_types,
        [ctypes.POINTER(ctypes.c_char), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmPsf.getAtomTypes()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_charges,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmPsf.getCharges()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_masses,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmPsf.getMasses()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_net_charge,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CharmmPsf.getNetCharge()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_total_mass,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CharmmPsf.getTotalMass()",
    )

    configure_status_function(
        lib().apo_charmm_psf_get_file_name,
        [ctypes.POINTER(ctypes.c_char), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmPsf.getFileName()",
    )

    _prototypes_initialized = True

    return


class CharmmPsf(_ApoObject):
    _destroy_function_name = "apo_charmm_psf_destroy"

    def __init__(self, path: FilePath) -> None:
        _initialize_prototypes()
        super().__init__()

        handle = ctypes.c_void_p()

        encoded_path: bytes = encode_path(path)
        c_path: ctypes.c_char_p = ctypes.c_char_p(encoded_path)

        lib().apo_charmm_psf_create(ctypes.byref(handle), c_path)

        if handle.value is None:
            raise RuntimeError(
                "apo_charmm_psf_create returned success but produced a NULL handle"
            )

        self._handle = handle

        return

    def getNumAtoms(self) -> int:
        _initialize_prototypes()

        c_num_atoms = ctypes.c_size_t()

        lib().apo_charmm_psf_get_num_atoms(ctypes.byref(c_num_atoms), self.handle)

        return int(c_num_atoms.value)

    def getNumBonds(self) -> int:
        _initialize_prototypes()

        c_num_bonds = ctypes.c_size_t()

        lib().apo_charmm_psf_get_num_bonds(ctypes.byref(c_num_bonds), self.handle)

        return int(c_num_bonds.value)

    def getNumAngles(self) -> int:
        _initialize_prototypes()

        c_num_angles = ctypes.c_size_t()

        lib().apo_charmm_psf_get_num_angles(ctypes.byref(c_num_angles), self.handle)

        return int(c_num_angles.value)

    def getNumDihedrals(self) -> int:
        _initialize_prototypes()

        c_num_dihedrals = ctypes.c_size_t()

        lib().apo_charmm_psf_get_num_dihedrals(
            ctypes.byref(c_num_dihedrals), self.handle
        )

        return int(c_num_dihedrals.value)

    def getNumImpropers(self) -> int:
        _initialize_prototypes()

        c_num_impropers = ctypes.c_size_t()

        lib().apo_charmm_psf_get_num_impropers(
            ctypes.byref(c_num_impropers), self.handle
        )

        return int(c_num_impropers.value)

    def getNumCrossTerms(self) -> int:
        _initialize_prototypes()

        c_num_cross_terms = ctypes.c_size_t()

        lib().apo_charmm_psf_get_num_cross_terms(
            ctypes.byref(c_num_cross_terms), self.handle
        )

        return int(c_num_cross_terms.value)

    def getSegmentIdentifiers(self) -> list[str]:
        _initialize_prototypes()

        num_segis: int = self.getNumAtoms()

        c_buffer_type = ctypes.c_char * (num_segis * 8)
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_segis * 8)

        lib().apo_charmm_psf_get_segment_identifiers(
            c_buffer, c_buffer_len, self.handle
        )

        raw_buffer = c_buffer.raw

        segis: list[str] = []
        for i in range(num_segis):
            start = i * 8
            stop = start + 8
            segi = raw_buffer[start:stop].decode("utf-8")
            segis.append(segi.strip())  # Remove whitespace

        return segis

    def getResidueIdentifiers(self) -> list[int]:
        _initialize_prototypes()

        num_resis = self.getNumAtoms()

        c_buffer_type = ctypes.c_int * num_resis
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_resis)

        lib().apo_charmm_psf_get_residue_identifiers(
            c_buffer, c_buffer_len, self.handle
        )

        resis: list[int] = []
        for i in range(num_resis):
            resis.append(int(c_buffer[i]))

        return resis

    def getResidueNames(self) -> list[str]:
        _initialize_prototypes()

        num_resns = self.getNumAtoms()

        c_buffer_type = ctypes.c_char * (num_resns * 8)
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_resns * 8)

        lib().apo_charmm_psf_get_residue_names(c_buffer, c_buffer_len, self.handle)

        raw_buffer = c_buffer.raw

        resns: list[str] = []
        for i in range(num_resns):
            start = i * 8
            stop = start + 8
            resn = raw_buffer[start:stop].decode("utf-8")
            resns.append(resn.strip())  # Remove whitespace

        return resns

    def getAtomNames(self) -> list[str]:
        _initialize_prototypes()

        num_names = self.getNumAtoms()

        c_buffer_type = ctypes.c_char * (num_names * 8)
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_names * 8)

        lib().apo_charmm_psf_get_atom_names(c_buffer, c_buffer_len, self.handle)

        raw_buffer = c_buffer.raw

        atom_names: list[str] = []
        for i in range(num_names):
            start = i * 8
            stop = start + 8
            atom_name = raw_buffer[start:stop].decode("utf-8")
            atom_names.append(atom_name.strip())  # Remove whitespace

        return atom_names

    def getAtomTypes(self) -> list[str]:
        _initialize_prototypes()

        num_types = self.getNumAtoms()

        c_buffer_type = ctypes.c_char * (num_types * 8)
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_types * 8)

        lib().apo_charmm_psf_get_atom_types(c_buffer, c_buffer_len, self.handle)

        raw_buffer = c_buffer.raw

        atom_types: list[str] = []
        for i in range(num_types):
            start = i * 8
            stop = start + 8
            atom_type = raw_buffer[start:stop].decode("utf-8")
            atom_types.append(atom_type.strip())  # Remove whitespace

        return atom_types

    def getCharges(self) -> list[float]:
        _initialize_prototypes()

        num_atoms = self.getNumAtoms()

        c_buffer_type = ctypes.c_double * num_atoms
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_atoms)

        lib().apo_charmm_psf_get_charges(c_buffer, c_buffer_len, self.handle)

        charges: list[float] = []
        for i in range(num_atoms):
            charges.append(float(c_buffer[i]))

        return charges

    def getMasses(self) -> list[float]:
        _initialize_prototypes()

        num_atoms = self.getNumAtoms()

        c_buffer_type = ctypes.c_double * num_atoms
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_atoms)

        lib().apo_charmm_psf_get_masses(c_buffer, c_buffer_len, self.handle)

        masses: list[float] = []
        for i in range(num_atoms):
            masses.append(float(c_buffer[i]))

        return masses

    def getNetCharge(self) -> float:
        _initialize_prototypes()

        c_net_charge = ctypes.c_double()

        lib().apo_charmm_psf_get_net_charge(ctypes.byref(c_net_charge), self.handle)

        return float(c_net_charge.value)

    def getTotalMass(self) -> float:
        _initialize_prototypes()

        c_total_mass = ctypes.c_double()

        lib().apo_charmm_psf_get_total_mass(ctypes.byref(c_total_mass), self.handle)

        return float(c_total_mass.value)

    def getFileName(self) -> str:
        _initialize_prototypes()

        # We are going to assume that the file name is less than 1024 chars
        c_buffer_type = ctypes.c_char * 1024
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(1024)

        lib().apo_charmm_psf_get_file_name(c_buffer, c_buffer_len, self.handle)

        file_name = c_buffer.raw.decode("utf-8")

        return file_name.rstrip(" ")
