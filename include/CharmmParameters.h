// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE
//

/* This class contains CHARMM parameters
 * It is filled up by reading the CHARMM .prm file
 *
 */
#pragma once

#include "CharmmPSF.h"

#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

struct BondKey {
public:
  BondKey(const std::string &iatom, const std::string &jatom)
      : atom1(iatom), atom2(jatom) {}

public:
  bool operator==(const BondKey &other) const {
    return this->atom1 == other.atom1 && this->atom2 == other.atom2;
  }

  friend bool operator<(const BondKey &first, const BondKey &second) {
    if (second.atom1 != first.atom1)
      return second.atom1 < first.atom1;
    else
      return second.atom2 < first.atom2;
  }

  friend std::ostream &operator<<(std::ostream &output, const BondKey &key) {
    output << key.atom1 << " " << key.atom2 << " ";
    return output;
  }

public:
  std::string atom1;
  std::string atom2;
};

class BondValues {
public:
  BondValues(void) : kb(0.0), b0(0.0) {}
  BondValues(const double kb, const double b0) : kb(kb), b0(b0) {}
  BondValues(const BondValues &other) = default;

public:
  friend std::ostream &operator<<(std::ostream &output, const BondValues &bv) {
    output << "(kb: " << bv.kb << ", b0:" << bv.b0 << ")\n";
    return output;
  }

public:
  double kb;
  double b0;
};

struct AngleKey {
public:
  AngleKey(const std::string &iatom, const std::string &jatom,
           const std::string &katom)
      : atom1(iatom), atom2(jatom), atom3(katom) {}

public:
  bool operator==(const AngleKey &other) const {
    return this->atom1 == other.atom1 && this->atom2 == other.atom2 &&
           this->atom3 == other.atom3;
  }

  friend bool operator<(const AngleKey &first, const AngleKey &second) {
    if (second.atom1 != first.atom1)
      return second.atom1 < first.atom1;
    else if (second.atom2 != first.atom2)
      return second.atom2 < first.atom2;
    else
      return second.atom3 < first.atom3;
  }

  friend std::ostream &operator<<(std::ostream &output, const AngleKey &key) {
    output << key.atom1 << " " << key.atom2 << " " << key.atom3 << " ";
    return output;
  }

public:
  std::string atom1;
  std::string atom2;
  std::string atom3;
};

struct AngleValues {
public:
  AngleValues(void) : kTheta(0.0), theta0(0.0) {}
  AngleValues(const double kTheta, const double theta0)
      : kTheta(kTheta), theta0(theta0) {}
  AngleValues(const AngleValues &other) = default;

public:
  friend std::ostream &operator<<(std::ostream &output, const AngleValues &av) {
    output << "(kTheta: " << av.kTheta << ", theta0: " << av.theta0 << ")\n";
    return output;
  }

public:
  double kTheta;
  double theta0;
};

struct DihedralKey {
public:
  DihedralKey(const std::string &iatom, const std::string &jatom,
              const std::string &katom, const std::string &latom)
      : atom1(iatom), atom2(jatom), atom3(katom), atom4(latom) {}

public:
  friend bool operator<(const DihedralKey &first, const DihedralKey &second) {
    if (second.atom1 != first.atom1)
      return second.atom1 < first.atom1;
    else if (second.atom2 != first.atom2)
      return second.atom2 < first.atom2;
    else if (second.atom3 != first.atom3)
      return second.atom3 < first.atom3;
    else
      return second.atom4 < first.atom4;
  }

  friend std::ostream &operator<<(std::ostream &output, const DihedralKey &dk) {
    output << dk.atom1 << " " << dk.atom2 << " " << dk.atom3 << " " << dk.atom4
           << " ";
    return output;
  }

  bool operator==(const DihedralKey &other) const {
    return this->atom1 == other.atom1 && this->atom2 == other.atom2 &&
           this->atom3 == other.atom3 && this->atom4 == other.atom4;
  }

public:
  std::string atom1;
  std::string atom2;
  std::string atom3;
  std::string atom4;
};

struct DihedralValues {
public:
  DihedralValues(void) : kChi(0.0), delta(0.0), n(0) {}
  DihedralValues(const double kChi, const int n, const double delta)
      : kChi(kChi), delta(delta), n(n) {}
  DihedralValues(const DihedralValues &other) = default;

public:
  friend std::ostream &operator<<(std::ostream &output,
                                  const DihedralValues &dv) {
    output << "(kChi: " << dv.kChi << ", delta : " << dv.delta
           << ", n : " << dv.n << ")\n";
    return output;
  }

public:
  double kChi;
  double delta;
  int n;
};

struct ImDihedralValues {
public:
  ImDihedralValues(void) : kPsi(0.0), psi0(0.0) {}
  ImDihedralValues(const double kPsi, const double psi0)
      : kPsi(kPsi), psi0(psi0) {}
  ImDihedralValues(const ImDihedralValues &other) = default;

public:
  friend std::ostream &operator<<(std::ostream &output,
                                  const ImDihedralValues &idv) {
    output << "(kPsi : " << idv.kPsi << ", psi0 : " << idv.psi0 << ")\n";
    return output;
  }

public:
  double kPsi;
  double psi0;
};

struct CmapKey {
public:
  CmapKey(const DihedralKey &d1, const DihedralKey &d2) : dih1(d1), dih2(d2) {}

public:
  DihedralKey dih1;
  DihedralKey dih2;
};

class VdwParameters {
public:
  VdwParameters(void) : epsilon(0.0), rmin_2(0.0) {}
  VdwParameters(const double eps, const double rm2)
      : epsilon(eps), rmin_2(rm2) {}

public:
  friend std::ostream &operator<<(std::ostream &output,
                                  const VdwParameters &vdw) {
    output << "(epsilon: " << vdw.epsilon << ", rmin_2: " << vdw.rmin_2
           << ")\n";
    return output;
  }

public:
  double epsilon;
  double rmin_2;
};

// JEG260812: This and CMAP need to be fixed
class NBFixParameters {
public:
  // NBFixParameters(std::string a1, std::string a2, double e, double r,
  //                 double e14, double r14)
  //     : atom1(a1), atom2(a2), emin(e), rmin(r), emin14(e14), rmin14(r14) {}
  // // copy constructor
  // NBFixParameters(const NBFixParameters &nbfix) = default;
  // NBFixParameters(const NBFixParameters &nbfix)
  //     : atom1(nbfix.atom1), atom2(nbfix.atom2), emin(nbfix.emin),
  //       rmin(nbfix.rmin), emin14(nbfix.emin14), rmin14(nbfix.rmin14) {}

  // // assignment operator
  // NBFixParameters &operator=(const NBFixParameters &nbfix) {
  //   atom1 = nbfix.atom1;
  //   atom2 = nbfix.atom2;
  //   emin = nbfix.emin;
  //   rmin = nbfix.rmin;
  //   emin14 = nbfix.emin14;
  //   rmin14 = nbfix.rmin14;
  //   return *this;
  // }

  friend std::ostream &operator<<(std::ostream &output,
                                  const NBFixParameters &nbfix) {
    output << "(" << nbfix.atom1 << "," << nbfix.atom2 << "," << nbfix.emin
           << "," << nbfix.rmin << "," << nbfix.emin14 << "," << nbfix.rmin14
           << ")\n";
    return output;
  }

public:
  std::string atom1, atom2;
  double emin, rmin, emin14, rmin14;
};

/** @brief Contains bonded interactions parameters and list
 * @todo improve doc
 *
 * Contains:
 * - paramsSize: vector of int, containing the number of types of bond,
 * Urey-Bradley, angle, dihedral and (in the future ?) CMAP
 * - paramsVal: vector of vectors, containing the parameters for each
 * interaction
 * - listsSize: vector of int, contains the number of interactions for each type
 * - listVal: vector of vectors of int, contains the atom indices for each
 * interaction
 */
struct BondedParamsAndLists {
public:
  BondedParamsAndLists(const std::vector<int> &pSize,
                       const std::vector<std::vector<float>> &pVal,
                       const std::vector<int> &lSize,
                       const std::vector<std::vector<int>> &lVal)
      : paramsSize(pSize), paramsVal(pVal), listsSize(lSize), listVal(lVal) {}

public:
  std::vector<int> paramsSize;
  std::vector<std::vector<float>> paramsVal;

  std::vector<int> listsSize;
  std::vector<std::vector<int>> listVal;
};

/** @brief Contains van der Waals interactions parameters and lists */
struct VdwParamsAndTypes {
public:
  VdwParamsAndTypes(const std::vector<float> &vdwParams,
                    const std::vector<float> &vdw14Params,
                    const std::vector<int> &vdwTypes,
                    const std::vector<int> &vdw14Types)
      : vdwParams(vdwParams), vdw14Params(vdw14Params), vdwTypes(vdwTypes),
        vdw14Types(vdw14Types) {}

public:
  std::vector<float> vdwParams, vdw14Params;
  std::vector<int> vdwTypes, vdw14Types;
};

/**
 * @brief Set of CHARMM parameters
 *
 * CHARMM parameters. Initialized from a filename (.prm, .str) or a list of
 * file names. Does not contain charges.
 */
class CharmmParameters {
public:
  CharmmParameters(void);
  /** @brief Constructor. Uses a single .prm or .str input file.
   * @param fileName String, .prm or .str file
   */
  CharmmParameters(const std::string &fileName);
  /** @brief Constructor. Uses a list of .prm and .str input files.
   * @param fileNames List of filename strings
   */
  CharmmParameters(const std::vector<std::string> &fileNames);

public: // Getters
  /** @brief Returns a map of the bond interactions. Key is a BondKey object
   * (couple of atom names), value is a BondValues object (kb and b0) */
  const std::map<BondKey, BondValues> &getBondParams(void) const;
  const std::map<AngleKey, BondValues> &getUreybParams(void) const;
  const std::map<AngleKey, AngleValues> &getAngleParams(void) const;
  const std::map<DihedralKey, std::vector<DihedralValues>> &
  getDihedralParams(void) const;
  const std::map<DihedralKey, ImDihedralValues> &getImproperParams(void) const;
  const std::map<std::string, VdwParameters> &getVdwParams(void) const;
  const std::map<std::string, VdwParameters> &getVdw14Params(void) const;

  /** @brief Get original names of prm files imported */
  const std::vector<std::string> &getPrmFileNames(void) const;

public:
  /** @brief Returns bonded interactions list and parameters as
   * BondedParamsAndLists object
   *
   */
  BondedParamsAndLists
  getBondedParamsAndLists(const std::shared_ptr<CharmmPSF> &psf) const;

  /** @brief Returns vdW interaction lists and parameters as VdwParamsAndTypes
   * object */
  VdwParamsAndTypes getVdwParamsAndTypes(std::shared_ptr<CharmmPSF> &psf) const;

  /**
   * @brief Read prm file
   *
   * @param[in] fileName Input file name
   */
  void readCharmmParameterFile(const std::string &fileName);

private:
  void parseBondRecord(const std::vector<std::string> &tokens,
                       const std::string &line, const std::string &fileName,
                       const std::size_t lineNumber);
  void parseAngleRecord(const std::vector<std::string> &tokens,
                        const std::string &line, const std::string &fileName,
                        const std::size_t lineNumber);
  void parseDihedralRecord(const std::vector<std::string> &tokens,
                           const std::string &line, const std::string &fileName,
                           const std::size_t lineNumber);
  void parseImproperRecord(const std::vector<std::string> &tokens,
                           const std::string &line, const std::string &fileName,
                           const std::size_t lineNumber);
  void parseNonbondedRecord(const std::vector<std::string> &tokens,
                            const std::string &line,
                            const std::string &fileName,
                            const std::size_t lineNumber);
  void parseNbfixRecord(const std::vector<std::string> &tokens,
                        const std::string &line, const std::string &fileName,
                        const std::size_t lineNumber);

private:
  std::map<BondKey, BondValues> m_BondParams;
  std::map<AngleKey, BondValues> m_UreybParams;
  std::map<AngleKey, AngleValues> m_AngleParams;
  std::map<DihedralKey, std::vector<DihedralValues>> m_DihedralParams;
  std::map<DihedralKey, ImDihedralValues> m_ImproperParams;

  // will be filled in the  vdw(14)Params
  std::map<std::tuple<std::string, std::string>, NBFixParameters> m_NbfixParams;

  std::map<std::string, VdwParameters> m_VdwParams;
  std::map<std::string, VdwParameters> m_Vdw14Params;

  /** @brief original prm file names */
  std::vector<std::string> m_PrmFileNames;
};
