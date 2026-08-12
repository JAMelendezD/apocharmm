// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#include "CharmmParameters.h"

#include "ApoCharmmError.h"
#include "str_utils.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

template <typename T>
T ParsePrmValue(const std::string &token, const std::string_view fieldName,
                const std::string_view recordType, const std::string &line,
                const std::string &fileName, const std::size_t lineNumber) {
  static_assert(std::is_same_v<T, double> || std::is_same_v<T, int>,
                "Unsupported CHARMM parameter value type");

  std::string normalizedToken = token;
  if constexpr (std::is_floating_point_v<T>) {
    std::replace(normalizedToken.begin(), normalizedToken.end(), 'D', 'E');
    std::replace(normalizedToken.begin(), normalizedToken.end(), 'd', 'E');
  }

  std::size_t parsedCharacters = 0;
  T value{};
  bool conversionSucceeded = false;

  try {
    if constexpr (std::is_same_v<T, double>)
      value = std::stod(normalizedToken, &parsedCharacters);
    else
      value = std::stoi(normalizedToken, &parsedCharacters);

    conversionSucceeded = true;
  } catch (const std::invalid_argument &) {
  } catch (const std::out_of_range &) {
  }

  bool isValid =
      conversionSucceeded && (parsedCharacters == normalizedToken.size());

  if constexpr (std::is_floating_point_v<T>)
    isValid = isValid && std::isfinite(value);

  APOCHARMM_REQUIRE(
      isValid, ApoCharmmErrorCode::Runtime,
      "Invalid " + std::string(fieldName) + " value \"" + token + "\" in " +
          std::string(recordType) + " parameter record in file \"" + fileName +
          "\" at line " + std::to_string(lineNumber) + ": " + line);

  return value;
}

std::string CleanPrmLine(const std::string &rawLine) {
  std::string line = rawLine;

  const std::size_t commentPosition = line.find('!');
  if (commentPosition != std::string::npos)
    line.erase(commentPosition);

  std::replace(line.begin(), line.end(), '\t', ' ');
  apo::trim_ip(line);
  apo::to_upper_ip(line);

  return line;
}

} // namespace

CharmmParameters::CharmmParameters(void)
    : m_BondParams(), m_UreybParams(), m_AngleParams(), m_DihedralParams(),
      m_ImproperParams(), m_NbfixParams(), m_VdwParams(), m_Vdw14Params(),
      m_PrmFileNames() {}

CharmmParameters::CharmmParameters(const std::string &fileName)
    : CharmmParameters() {
  m_PrmFileNames.push_back(fileName);
  this->readCharmmParameterFile(fileName);
}

CharmmParameters::CharmmParameters(const std::vector<std::string> &fileNames)
    : CharmmParameters() {
  for (const std::string &fileName : fileNames) {
    m_PrmFileNames.push_back(fileName);
    this->readCharmmParameterFile(fileName);
  }
}

const std::map<BondKey, BondValues> &
CharmmParameters::getBondParams(void) const {
  return m_BondParams;
}

const std::map<AngleKey, BondValues> &
CharmmParameters::getUreybParams(void) const {
  return m_UreybParams;
}

const std::map<AngleKey, AngleValues> &
CharmmParameters::getAngleParams(void) const {
  return m_AngleParams;
}

const std::map<DihedralKey, std::vector<DihedralValues>> &
CharmmParameters::getDihedralParams(void) const {
  return m_DihedralParams;
}

const std::map<DihedralKey, ImDihedralValues> &
CharmmParameters::getImproperParams(void) const {
  return m_ImproperParams;
}

const std::map<std::string, VdwParameters> &
CharmmParameters::getVdwParams(void) const {
  return m_VdwParams;
}

const std::map<std::string, VdwParameters> &
CharmmParameters::getVdw14Params(void) const {
  return m_Vdw14Params;
}

const std::vector<std::string> &CharmmParameters::getPrmFileNames(void) const {
  return m_PrmFileNames;
}

BondedParamsAndLists CharmmParameters::getBondedParamsAndLists(
    const std::shared_ptr<CharmmPSF> &psf) const {
  std::vector<int> paramsSize;
  std::vector<std::vector<float>> paramsVal;

  std::vector<int> listsSize;
  std::vector<std::vector<int>> listVal;

  const std::vector<std::string> &atomTypes = psf->getAtomTypes();
  const std::vector<std::string> &atomNames = psf->getAtomNames();
  const std::vector<Bond> &bonds = psf->getBonds();
  const std::vector<Angle> &angles = psf->getAngles();
  const std::vector<Dihedral> &dihedrals = psf->getDihedrals();
  const std::vector<Dihedral> &impropers = psf->getImpropers();
  const std::vector<CrossTerm> &cmaps = psf->getCrossTerms();

  std::vector<BondKey> bondKeysPresent;
  std::vector<AngleKey> ureybKeysPresent;
  std::vector<AngleKey> angleKeysPresent;
  std::vector<DihedralKey> dihedralKeysPresent;
  std::vector<DihedralKey> improperKeysPresent;

  for (int bond = 0; bond < psf->getNumBonds(); bond++) {
    std::string atom1 = atomTypes[bonds[bond].iatom];
    std::string atom2 = atomTypes[bonds[bond].jatom];
    if (atom1 > atom2)
      std::swap(atom1, atom2);
    auto key = BondKey(atom1, atom2);

    if (m_BondParams.count(key)) {
      auto findResult =
          std::find(bondKeysPresent.begin(), bondKeysPresent.end(), key);
      if (findResult == std::end(bondKeysPresent)) {
        bondKeysPresent.push_back(key);
        const BondValues &value = m_BondParams.at(key);
        paramsVal.push_back(
            {static_cast<float>(value.b0), static_cast<float>(value.kb)});
      }
      findResult =
          std::find(bondKeysPresent.begin(), bondKeysPresent.end(), key);
      int bondType = findResult - std::begin(bondKeysPresent);
      listVal.push_back({bonds[bond].iatom, bonds[bond].jatom, bondType, 13});
    } else {
      std::stringstream tmpexc;
      tmpexc << "bond not found " << bond << " " << key << " "
             << bonds[bond].iatom << " " << bonds[bond].jatom << "\n";
      throw std::invalid_argument(tmpexc.str());
    }
  }
  paramsSize.push_back(bondKeysPresent.size());
  listsSize.push_back(listVal.size());

  // Most angles do not have a Urey-Bradley contribution, so only register terms
  // whose force constant is nonzero.
  for (int angle = 0; angle < psf->getNumAngles(); angle++) {
    std::string atom1 = atomTypes[angles[angle].iatom];
    std::string atom2 = atomTypes[angles[angle].jatom];
    std::string atom3 = atomTypes[angles[angle].katom];
    if (atom1 > atom3)
      std::swap(atom1, atom3);
    auto key = AngleKey(atom1, atom2, atom3);

    if (m_UreybParams.count(key)) {
      const BondValues &value = m_UreybParams.at(key);
      if (std::abs(value.kb) <= 0.01)
        continue;

      auto findResult =
          std::find(ureybKeysPresent.begin(), ureybKeysPresent.end(), key);
      if (findResult == ureybKeysPresent.end()) {
        ureybKeysPresent.push_back(key);
        paramsVal.push_back(
            {static_cast<float>(value.b0), static_cast<float>(value.kb)});
      }
      findResult =
          std::find(ureybKeysPresent.begin(), ureybKeysPresent.end(), key);
      int ureybType = findResult - ureybKeysPresent.begin();
      listVal.push_back(
          {angles[angle].iatom, angles[angle].katom, ureybType, 13});
    } else {
      std::stringstream tmpexc;
      tmpexc << "Ureyb not found " << angle << " " << key << "\n";
      throw std::invalid_argument(tmpexc.str());
    }
  }
  paramsSize.push_back(ureybKeysPresent.size());
  listsSize.push_back(listVal.size() - listsSize[0]);

  for (int angle = 0; angle < psf->getNumAngles(); angle++) {
    std::string atom1 = atomTypes[angles[angle].iatom];
    std::string atom2 = atomTypes[angles[angle].jatom];
    std::string atom3 = atomTypes[angles[angle].katom];
    if (atom1 > atom3)
      std::swap(atom1, atom3);
    auto key = AngleKey(atom1, atom2, atom3);

    if (m_AngleParams.count(key)) {
      auto findResult =
          std::find(angleKeysPresent.begin(), angleKeysPresent.end(), key);
      if (findResult == angleKeysPresent.end()) {
        angleKeysPresent.push_back(key);
        const AngleValues &value = m_AngleParams.at(key);
        paramsVal.push_back({static_cast<float>(value.theta0),
                             static_cast<float>(value.kTheta)});
      }
      findResult =
          std::find(angleKeysPresent.begin(), angleKeysPresent.end(), key);
      int angleType = findResult - angleKeysPresent.begin();
      listVal.push_back({angles[angle].iatom, angles[angle].jatom,
                         angles[angle].katom, angleType, 13, 13});

    } else {
      std::stringstream tmpexc;
      tmpexc << "Angle not found " << angle << " " << key << "\n";
      throw std::invalid_argument(tmpexc.str());
    }
  }
  paramsSize.push_back(angleKeysPresent.size());
  listsSize.push_back(listVal.size() - listsSize[0] - listsSize[1]);

  int dihedralParamsPresent = 0;
  int startParamDihedral = paramsVal.size();
  std::map<DihedralKey, int> indexOfKeyInParamsVal;
  const double pi_180 = std::acos(-1) / 180.0;
  for (int dihedral = 0; dihedral < psf->getNumDihedrals(); dihedral++) {
    std::string atom1 = atomTypes[dihedrals[dihedral].iatom];
    std::string atom2 = atomTypes[dihedrals[dihedral].jatom];
    std::string atom3 = atomTypes[dihedrals[dihedral].katom];
    std::string atom4 = atomTypes[dihedrals[dihedral].latom];

    if (atom1 > atom4) {
      std::swap(atom1, atom4);
      std::swap(atom2, atom3);
    }
    if ((atom1 == atom4) && (atom2 > atom3))
      std::swap(atom2, atom3);
    auto key = DihedralKey(atom1, atom2, atom3, atom4);

    if (m_DihedralParams.count(key)) {
      auto findResult = std::find(dihedralKeysPresent.begin(),
                                  dihedralKeysPresent.end(), key);
      if (findResult == dihedralKeysPresent.end()) {
        const std::vector<DihedralValues> &values = m_DihedralParams.at(key);

        for (std::size_t i = 0; i < values.size(); i++) {
          if (i > 0)
            paramsVal[paramsVal.size() - 1][0] *= -1;
          else {
            indexOfKeyInParamsVal[key] =
                static_cast<int>(paramsVal.size() - startParamDihedral);
          }

          const DihedralValues &value = values[i];
          const double cpsin = std::sin(value.delta * pi_180);
          const double cpcos = std::cos(value.delta * pi_180);
          paramsVal.push_back(
              {static_cast<float>(value.n), static_cast<float>(value.kChi),
               static_cast<float>(cpsin), static_cast<float>(cpcos)});
          dihedralParamsPresent++;
        }
        dihedralKeysPresent.push_back(key);
      }
      int dihedralType = indexOfKeyInParamsVal[key];
      listVal.push_back({dihedrals[dihedral].iatom, dihedrals[dihedral].jatom,
                         dihedrals[dihedral].katom, dihedrals[dihedral].latom,
                         dihedralType, 13, 13, 13});
    } else {
      if (atom2 > atom3)
        std::swap(atom2, atom3);
      key = DihedralKey("X", atom2, atom3, "X");

      if (m_DihedralParams.count(key)) {
        auto findResult = std::find(dihedralKeysPresent.begin(),
                                    dihedralKeysPresent.end(), key);
        if (findResult == dihedralKeysPresent.end()) {
          const std::vector<DihedralValues> &values = m_DihedralParams.at(key);

          for (std::size_t i = 0; i < values.size(); i++) {
            if (i > 0)
              paramsVal[paramsVal.size() - 1][0] *= -1;
            else {
              indexOfKeyInParamsVal[key] =
                  static_cast<int>(paramsVal.size() - startParamDihedral);
            }

            const DihedralValues &value = values[i];
            const double cpsin = std::sin(value.delta * pi_180);
            const double cpcos = std::cos(value.delta * pi_180);
            paramsVal.push_back(
                {static_cast<float>(value.n), static_cast<float>(value.kChi),
                 static_cast<float>(cpsin), static_cast<float>(cpcos)});
            dihedralParamsPresent++;
          }
          dihedralKeysPresent.push_back(key);
        }
        int dihedralType = indexOfKeyInParamsVal[key];
        listVal.push_back({dihedrals[dihedral].iatom, dihedrals[dihedral].jatom,
                           dihedrals[dihedral].katom, dihedrals[dihedral].latom,
                           dihedralType, 13, 13, 13});
      } else {
        std::stringstream tmpexc;
        tmpexc << "dihedral not found " << dihedral << " "
               << atomTypes[dihedrals[dihedral].iatom] << " "
               << atomTypes[dihedrals[dihedral].jatom] << " "
               << atomTypes[dihedrals[dihedral].katom] << " "
               << atomTypes[dihedrals[dihedral].latom] << "\t"
               << atomNames[dihedrals[dihedral].iatom] << " "
               << atomNames[dihedrals[dihedral].jatom] << " "
               << atomNames[dihedrals[dihedral].katom] << " "
               << atomNames[dihedrals[dihedral].latom] << "\n";
        throw std::invalid_argument(tmpexc.str());
      }
    }
  }
  paramsSize.push_back(dihedralParamsPresent);
  listsSize.push_back(listVal.size() - listsSize[0] - listsSize[1] -
                      listsSize[2]);

  for (int improper = 0; improper < psf->getNumImpropers(); improper++) {
    std::string atom1 = atomTypes[impropers[improper].iatom];
    std::string atom2 = atomTypes[impropers[improper].jatom];
    std::string atom3 = atomTypes[impropers[improper].katom];
    std::string atom4 = atomTypes[impropers[improper].latom];
    if (atom1 > atom4) {
      std::swap(atom1, atom4);
      std::swap(atom2, atom3);
    }
    if ((atom1 == atom4) && (atom2 > atom3))
      std::swap(atom2, atom3);
    auto key = DihedralKey(atom1, atom2, atom3, atom4);

    if (m_ImproperParams.count(key)) {
      auto findResult = std::find(improperKeysPresent.begin(),
                                  improperKeysPresent.end(), key);
      if (findResult == improperKeysPresent.end()) {
        improperKeysPresent.push_back(key);
        const ImDihedralValues &value = m_ImproperParams.at(key);
        paramsVal.push_back({static_cast<float>(value.psi0),
                             static_cast<float>(value.kPsi), 0.0f, 1.0f});
      }
      findResult = std::find(improperKeysPresent.begin(),
                             improperKeysPresent.end(), key);
      int improperType = findResult - improperKeysPresent.begin();
      listVal.push_back({impropers[improper].iatom, impropers[improper].jatom,
                         impropers[improper].katom, impropers[improper].latom,
                         improperType, 13, 13, 13});
    } else {
      key = DihedralKey(atom1, "X", "X", atom4);

      if (m_ImproperParams.count(key)) {
        auto findResult = std::find(improperKeysPresent.begin(),
                                    improperKeysPresent.end(), key);
        if (findResult == improperKeysPresent.end()) {
          improperKeysPresent.push_back(key);
          const ImDihedralValues &value = m_ImproperParams.at(key);
          paramsVal.push_back({static_cast<float>(value.psi0),
                               static_cast<float>(value.kPsi), 0.0f, 1.0f});
        }
        findResult = std::find(improperKeysPresent.begin(),
                               improperKeysPresent.end(), key);
        int improperType = findResult - improperKeysPresent.begin();
        listVal.push_back({impropers[improper].iatom, impropers[improper].jatom,
                           impropers[improper].katom, impropers[improper].latom,
                           improperType, 13, 13, 13});
      } else {
        std::stringstream tmpexc;
        tmpexc << "improper not found " << improper << " "
               << atomTypes[impropers[improper].iatom] << " "
               << atomTypes[impropers[improper].jatom] << " "
               << atomTypes[impropers[improper].katom] << " "
               << atomTypes[impropers[improper].latom] << "\t"
               << atomNames[impropers[improper].iatom] << " "
               << atomNames[impropers[improper].jatom] << " "
               << atomNames[impropers[improper].katom] << " "
               << atomNames[impropers[improper].latom] << "\n";
        throw std::invalid_argument(tmpexc.str());
      }
    }
  }
  paramsSize.push_back(improperKeysPresent.size());
  listsSize.push_back(listVal.size() - listsSize[0] - listsSize[1] -
                      listsSize[2] - listsSize[3]);

  for (int i = 0; i < psf->getNumCrossTerms(); ++i) {
    auto cmap = cmaps[i];
    // std::cout << cmap.atomi1 << " " << cmap.atomj1 << " " << cmap.atomk1 <<
    // "
    // "
    //           << cmap.atoml1 << " " << cmap.atomi2 << " " << cmap.atomj2 <<
    //           "
    //           "
    //           << cmap.atomk2 << " " << cmap.atoml2 << "\n";
    // std::cout << atomTypes[cmap.atomi1] << " " << atomTypes[cmap.atomj1] <<
    // "
    // "
    //           << atomTypes[cmap.atomk1] << " " << atomTypes[cmap.atoml1] <<
    //           "
    //           "
    //           << atomTypes[cmap.atomi2] << " " << atomTypes[cmap.atomj2] <<
    //           "
    //           "
    //           << atomTypes[cmap.atomk2] << " " << atomTypes[cmap.atoml2]
    //           << "\n";
    auto dihe1 = DihedralKey(atomTypes[cmap.iatom1], atomTypes[cmap.jatom1],
                             atomTypes[cmap.katom1], atomTypes[cmap.latom1]);
    auto dihe2 = DihedralKey(atomTypes[cmap.iatom2], atomTypes[cmap.jatom2],
                             atomTypes[cmap.katom2], atomTypes[cmap.latom2]);
    // std::cout << dihe1 << " " << dihe2 << "\n";
    auto key = CmapKey(dihe1, dihe2);
    // auto key = CmapKey(atomTypes[cmap.atom1], atomTypes[cmap.atom2],
    //                    atomTypes[cmap.atom3], atomTypes[cmap.atom4],
    //                    atomTypes[cmap.atom5]);
    // if (cmapParams.count(key)) {
    //   auto findResult =
    //       std::find(cmapKeysPresent.begin(), cmapKeysPresent.end(), key);
    //   if (findResult == cmapKeysPresent.end()) {
    //     cmapKeysPresent.push_back(key);
    //     auto value = cmapParams[key];
    //     paramsVal.push_back(value);
    //   }
    //   findResult =
    //       std::find(cmapKeysPresent.begin(), cmapKeysPresent.end(), key);
    //   int cmapType = findResult - cmapKeysPresent.begin();
    //   listVal.push_back({cmap.atom1, cmap.atom2, cmap.atom3, cmap.atom4,
    //                      cmap.atom5, cmapType, 13, 13, 13, 13});
    // } else {
    //   std::stringstream tmpexc;
    //   tmpexc << "cmap not found " << i << " " << key << "\n";
    //   throw std::invalid_argument(tmpexc.str());
    // }
  }

  // CMAP are currently not being used
  paramsSize.push_back(0);
  listsSize.push_back(0);
  return BondedParamsAndLists(paramsSize, paramsVal, listsSize, listVal);
}

VdwParamsAndTypes
CharmmParameters::getVdwParamsAndTypes(std::shared_ptr<CharmmPSF> &psf) const {
  std::vector<float> psfVdwParams, psfVdw14Params;
  std::vector<int> psfVdwTypes, psfVdw14Types;
  std::set<std::string> vdwAtomTypesMap, vdw14AtomTypesMap;

  for (const std::string &atomType : psf->getAtomTypes()) {
    vdwAtomTypesMap.insert(atomType);
    auto findResult = m_Vdw14Params.find(atomType);
    if (findResult != m_Vdw14Params.end())
      vdw14AtomTypesMap.insert(atomType);
  }

  std::vector<std::string> vdwAtomTypes(vdwAtomTypesMap.begin(),
                                        vdwAtomTypesMap.end());
  std::vector<std::string> vdw14AtomTypes(vdw14AtomTypesMap.begin(),
                                          vdw14AtomTypesMap.end());

  for (std::size_t i = 0; i < vdwAtomTypes.size(); i++) {
    for (std::size_t j = 0; j <= i; j++) {
      std::string iType = vdwAtomTypes[i];
      std::string jType = vdwAtomTypes[j];

      double epsilon, rmin;

      std::tuple<std::string, std::string> nbfixKey{jType, iType};
      if (m_NbfixParams.find(nbfixKey) != m_NbfixParams.end()) {
        const NBFixParameters &nbfix = m_NbfixParams.at(nbfixKey);
        epsilon = nbfix.emin;
        rmin = nbfix.rmin;
      } else {
        const VdwParameters &iParameters = m_VdwParams.at(iType);
        const VdwParameters &jParameters = m_VdwParams.at(jType);

        const double epsilonI = iParameters.epsilon;
        const double epsilonJ = jParameters.epsilon;
        const double rmin_2I = iParameters.rmin_2;
        const double rmin_2J = jParameters.rmin_2;

        epsilon = std::sqrt(epsilonI * epsilonJ);
        rmin = rmin_2I + rmin_2J;
      }

      double c12 = epsilon * std::pow(rmin, 12);
      double c6 = 2 * epsilon * std::pow(rmin, 6);

      psfVdwParams.push_back(static_cast<float>(c6));
      psfVdwParams.push_back(static_cast<float>(c12));
    }
  }

  for (std::size_t i = 0; i < vdwAtomTypes.size(); i++) {
    for (std::size_t j = 0; j <= i; j++) {
      std::string iType = vdwAtomTypes[i];
      std::string jType = vdwAtomTypes[j];

      double epsilon, rmin;

      std::tuple<std::string, std::string> nbfixKey{jType, iType};
      if (m_NbfixParams.find(nbfixKey) != m_NbfixParams.end()) {
        const NBFixParameters &nbfix = m_NbfixParams.at(nbfixKey);
        epsilon = nbfix.emin;
        rmin = nbfix.rmin;
      } else {
        const VdwParameters &iParameters = m_VdwParams.at(iType);
        const VdwParameters &jParameters = m_VdwParams.at(jType);

        double epsilonI = iParameters.epsilon;
        double epsilonJ = jParameters.epsilon;
        double rmin_2I = iParameters.rmin_2;
        double rmin_2J = jParameters.rmin_2;

        if (std::find(vdw14AtomTypes.begin(), vdw14AtomTypes.end(), iType) !=
            vdw14AtomTypes.end()) {
          const VdwParameters &i14Parameters = m_Vdw14Params.at(iType);
          epsilonI = i14Parameters.epsilon;
          rmin_2I = i14Parameters.rmin_2;
        }

        if (std::find(vdw14AtomTypes.begin(), vdw14AtomTypes.end(), jType) !=
            vdw14AtomTypes.end()) {
          const VdwParameters &j14Parameters = m_Vdw14Params.at(jType);
          epsilonJ = j14Parameters.epsilon;
          rmin_2J = j14Parameters.rmin_2;
        }

        epsilon = std::sqrt(epsilonI * epsilonJ);
        rmin = rmin_2I + rmin_2J;
      }

      double c12 = epsilon * std::pow(rmin, 12);
      double c6 = 2 * epsilon * std::pow(rmin, 6);

      psfVdw14Params.push_back(static_cast<float>(c6));
      psfVdw14Params.push_back(static_cast<float>(c12));
    }
  }

  int index = 0;
  for (const std::string &atomType : psf->getAtomTypes()) {
    auto result = std::find(vdwAtomTypes.begin(), vdwAtomTypes.end(), atomType);
    int pos = result - vdwAtomTypes.begin();
    psfVdwTypes.push_back(pos);
    psfVdw14Types.push_back(pos);
    index++;
  }

  return VdwParamsAndTypes(psfVdwParams, psfVdw14Params, psfVdwTypes,
                           psfVdw14Types);
}

void CharmmParameters::readCharmmParameterFile(const std::string &fileName) {
  enum class Section {
    NONE,
    ATOMS,
    BONDS,
    ANGLES,
    DIHEDRALS,
    IMPROPERS,
    CMAP,
    NONBONDED,
    NBFIX,
    HBOND
  };

  std::ifstream prmFile(fileName);
  APOCHARMM_REQUIRE(prmFile.is_open(), ApoCharmmErrorCode::Runtime,
                    "Failed to open CHARMM parameter file \"" + fileName +
                        "\"");

  std::size_t lineNumber = 0;
  std::string line = "";

  const std::string upperFileName = apo::to_upper(fileName);
  if (upperFileName.find("TOPPAR") != std::string::npos) {
    bool parameterBlockFound = false;

    while (std::getline(prmFile, line)) {
      lineNumber++;
      line = CleanPrmLine(line);

      if (line.empty() || (line.front() == '*'))
        continue;

      if ((line.find("READ") != std::string::npos) &&
          (line.find("PARA") != std::string::npos)) {
        parameterBlockFound = true;
        break;
      }
    }

    APOCHARMM_REQUIRE(parameterBlockFound, ApoCharmmErrorCode::Runtime,
                      "CHARMM parameter block was not found in file \"" +
                          fileName + "\"");
  }

  Section section = Section::NONE;

  while (std::getline(prmFile, line)) {
    lineNumber++;
    line = CleanPrmLine(line);

    if (line.empty() || (line.front() == '*'))
      continue;

    std::vector<std::string> tokens = apo::split(line);
    const std::string &keyword = tokens.front();

    if (keyword == "ATOMS") {
      section = Section::ATOMS;
      continue;
    }

    if (keyword == "BONDS") {
      section = Section::BONDS;
      continue;
    }

    if (keyword == "ANGLES") {
      section = Section::ANGLES;
      continue;
    }

    if (keyword == "DIHEDRALS") {
      section = Section::DIHEDRALS;
      continue;
    }

    if (keyword.rfind("IMPR", 0) == 0) {
      section = Section::IMPROPERS;
      continue;
    }

    if (keyword == "CMAP") {
      section = Section::CMAP;
      continue;
    }

    if (keyword == "NONBONDED") {
      section = Section::NONBONDED;
      const std::size_t headerLineNumber = lineNumber;

      while (!tokens.empty() && (tokens.back() == "-")) {
        bool continuationFound = false;

        while (std::getline(prmFile, line)) {
          lineNumber++;
          line = CleanPrmLine(line);

          if (line.empty() || (line.front() == '*'))
            continue;

          tokens = apo::split(line);
          continuationFound = true;
          break;
        }

        APOCHARMM_REQUIRE(continuationFound, ApoCharmmErrorCode::Runtime,
                          "Unexpected end of file while reading the NONBONDED "
                          "header in file \"" +
                              fileName + "\" beginning at line " +
                              std::to_string(headerLineNumber));
      }

      continue;
    }

    if (keyword == "NBFIX") {
      section = Section::NBFIX;
      continue;
    }

    if (keyword == "HBOND") {
      section = Section::HBOND;
      continue;
    }

    if (keyword == "END") {
      section = Section::NONE;
      continue;
    }

    switch (section) {
    case Section::BONDS:
      this->parseBondRecord(tokens, line, fileName, lineNumber);
      break;
    case Section::ANGLES:
      this->parseAngleRecord(tokens, line, fileName, lineNumber);
      break;
    case Section::DIHEDRALS:
      this->parseDihedralRecord(tokens, line, fileName, lineNumber);
      break;
    case Section::IMPROPERS:
      this->parseImproperRecord(tokens, line, fileName, lineNumber);
      break;
    case Section::NONBONDED:
      this->parseNonbondedRecord(tokens, line, fileName, lineNumber);
      break;
    case Section::NBFIX:
      this->parseNbfixRecord(tokens, line, fileName, lineNumber);
      break;
    case Section::NONE:
    case Section::ATOMS:
    case Section::CMAP:
    case Section::HBOND:
      break;
    }
  }

  APOCHARMM_REQUIRE(!prmFile.bad(), ApoCharmmErrorCode::Runtime,
                    "Failed while reading CHARMM parameter file \"" + fileName +
                        "\"");

  return;
}

/*
void CharmmParameters::readCharmmParameterFile(const std::string &fileName) {
  enum class State {
    NONE,
    ATOMS,
    BONDS,
    ANGLES,
    DIHEDRALS,
    IMPROPERS,
    CMAP,
    NONBONDED,
    NBFIX,
    HBOND
  };
  State state = State::NONE;
  std::vector<std::string> tokens;

  enum class FileType { NONE, PAR, TOPPAR };
  FileType fileType = FileType::NONE;
  std::size_t pos = fileName.find_last_of('/');
  if (pos != std::string::npos) {
    std::size_t topparPos = fileName.find("toppar", pos);
    if (topparPos != std::string::npos)
      fileType = FileType::TOPPAR;
    else
      fileType = FileType::PAR;
  }

  std::ifstream prmFile(fileName);
  std::string line;
  if (!prmFile.is_open()) {
    // std::cerr << "ERROR: Cannot open the file " << fileName << "\nExiting\n";
    throw std::invalid_argument("ERROR: Cannot open the file " + fileName +
                                "\nExiting\n");
    exit(0);
  }

  // If the file is toppar, skip the rtf portion
  // if (fileName.find_first_of("toppar") == 0 && fileName.find_last_of(".str")
  // != std::string::npos ) {

  if (fileType == FileType::TOPPAR) {

    while (!prmFile.eof()) {
      std::getline(prmFile, line);
      // std::cout << "toppar --" << line << "\n";
      apo::trim_ip(line);
      std::size_t pos = line.find_first_of('!');
      line = line.substr(0, pos);
      apo::trim_ip(line);
      apo::to_upper_ip(line);
      if (line.find_first_of('*') == 0 || line.find_first_of('!') == 0 ||
          line.size() == 0) {
        // Skip the line
      } else {
        if (line.find("READ") != std::string::npos &&
            line.find("PARA") != std::string::npos) {
          break;
        }
      }
    }
  }

  const float pi_180 = std::acos(-1) / 180.0;
  while (!prmFile.eof()) {
    std::getline(prmFile, line);
    apo::trim_ip(line);
    std::size_t pos = line.find_first_of('!');
    line = line.substr(0, pos);
    apo::trim_ip(line);
    apo::to_upper_ip(line);

    // std::cout << line << "\n";
    if (line.find_first_of('*') == 0 || line.find_first_of('!') == 0 ||
        line.size() == 0) {
      // Skip the line
    } else {
      if (line.find("ATOMS") == 0)
        state = State::ATOMS;
      if (line.find("BONDS") == 0)
        state = State::BONDS;
      if (line.find("ANGLES") == 0)
        state = State::ANGLES;
      if (line.find("DIHEDRALS") == 0)
        state = State::DIHEDRALS;
      if (line.find("IMPR") == 0)
        state = State::IMPROPERS;
      if (line.find("CMAP") == 0)
        state = State::CMAP;
      if (line.find("NONBONDED") == 0)
        state = State::NONBONDED;
      if (line.find("END") == 0)
        state = State::NONE;
      if (line.find("HBOND") == 0)
        state = State::HBOND;
      if (line.find("NBFIX") == 0)
        state = State::NBFIX;

      if (state == State::BONDS) {
        tokens = apo::split(line);
        if (tokens.size() >= 4) {
          if (tokens[0] > tokens[1])
            std::swap(tokens[0], tokens[1]);
          bondParams.insert(
              {BondKey(tokens[0], tokens[1]),
               BondValues(std::stof(tokens[2]), std::stof(tokens[3]))});
        }
      }
      if (state == State::ANGLES) {
        tokens = apo::split(line);
        // std::cout << "Tokens size : " << tokens.size() << "\n";
        if (tokens.size() >= 5) {
          if (tokens[0] > tokens[2])
            std::swap(tokens[0], tokens[2]);
          angleParams.insert({AngleKey(tokens[0], tokens[1], tokens[2]),
                              AngleValues(std::stof(tokens[3]),
                                          pi_180 * std::stof(tokens[4]))});
          // try {
          //   // std::cout << tokens[5] << "\n";
          //   float kub = std::stof(tokens[5]);
          //   float s0 = std::stof(tokens[6]);
          //   ureybParams.insert({AngleKey(tokens[0], tokens[1], tokens[2]),
          //                       BondValues(kub, s0)});
          //   // if (tokens[0] == "CT1" && tokens[2]=="CT2") std::cout << line
<<
          //   // "\n";
          // } catch (const std::exception &e) {
          //   ureybParams.insert({AngleKey(tokens[0], tokens[1], tokens[2]),
          //                       BondValues(0.0, 0.0)});
          //   // std::cout << line << " "  << e.what() << "\n";
          // }
          //
          if (tokens.size() > 5) {
            // std::cout << tokens[5] << "\n";
            float kub = std::stof(tokens[5]);
            float s0 = std::stof(tokens[6]);
            ureybParams.insert({AngleKey(tokens[0], tokens[1], tokens[2]),
                                BondValues(kub, s0)});
          } else {
            ureybParams.insert({AngleKey(tokens[0], tokens[1], tokens[2]),
                                BondValues(0.0, 0.0)});
          }

          //
        }
      }

      if (state == State::DIHEDRALS) {
        tokens = apo::split(line);
        if (tokens.size() >= 7) {
          if (tokens[0] > tokens[3]) {
            std::swap(tokens[0], tokens[3]);
            std::swap(tokens[2], tokens[1]);
          } else if ((tokens[0] == tokens[3]) && (tokens[1] > tokens[2]))
            std::swap(tokens[2], tokens[1]);

          auto key = DihedralKey(tokens[0], tokens[1], tokens[2], tokens[3]);
          auto elem = DihedralValues(std::stof(tokens[4]), std::stoi(tokens[5]),
                                     std::stof(tokens[6]));
          dihedralParams[key].push_back(elem);
          // std::cout << "new size " << dihedralParams[key].size() << "\n";
        }
      }

      if (state == State::IMPROPERS) {
        tokens = apo::split(line);
        if (tokens.size() >= 7) {
          if (tokens[0] > tokens[3]) {
            std::swap(tokens[0], tokens[3]);
            std::swap(tokens[2], tokens[1]);
          } else if ((tokens[0] == tokens[3]) && (tokens[1] > tokens[2]))
            std::swap(tokens[2], tokens[1]);
          improperParams.insert(
              {DihedralKey(tokens[0], tokens[1], tokens[2], tokens[3]),
               ImDihedralValues(std::stof(tokens[4]), std::stof(tokens[6]))});
        }
      }

      if (state == State::NONBONDED) {
        tokens = apo::split(line);
        if (tokens[0] == "NONBONDED") {
          while (*(tokens.end() - 1) == "-") {
            std::getline(prmFile, line);
            apo::ltrim_ip(line);
            std::vector<std::string> tokens1 = apo::split(line);
            tokens.insert(tokens.end(), tokens1.begin(), tokens1.end());
          }
        } else {
          apo::trim_ip(line);
          std::size_t pos = line.find_first_of('!');
          line = line.substr(0, pos);
          apo::trim_ip(line);
          tokens = apo::split(line);
          // if (tokens[0] == "HGA2")
          //   std::cout << line << " " << tokens.size() << "\n";
          //  std::cout << line << " " << tokens.size() << "\n";
          if (tokens.size() == 4) {
            if (vdwParams.find(tokens[0]) == vdwParams.end()) {
              vdwParams.insert(
                  {tokens[0],
                   VdwParameters(std::stod(tokens[2]), std::stod(tokens[3]))});
            } else {
              std::cerr << "Duplicate entry for " << tokens[0] << "\n";
              vdwParams[tokens[0]] =
                  VdwParameters(std::stod(tokens[2]), std::stod(tokens[3]));
              // throw std::invalid_argument("Duplicate entry for " + tokens[0]
              // +
              //                             "\n");
              // exit(0);
            }
            // vdwParams.insert({tokens[0], VdwParameters(std::stod(tokens[2]),
            // std::stod(tokens[3]))}); } else if (tokens.size() == 7) { if
(vdwParams.find(tokens[0]) == vdwParams.end()) { vdwParams.insert( {tokens[0],
                   VdwParameters(std::stod(tokens[2]), std::stod(tokens[3]))});
              vdw14Params.insert(
                  {tokens[0],
                   VdwParameters(std::stod(tokens[5]), std::stod(tokens[6]))});
            } else {
              std::cerr << "Duplicate entry for " << tokens[0] << "\n";
              vdwParams[tokens[0]] =
                  VdwParameters(std::stod(tokens[2]), std::stod(tokens[3]));
              vdw14Params[tokens[0]] =
                  VdwParameters(std::stod(tokens[5]), std::stod(tokens[6]));
              // throw std::invalid_argument("Duplicate entry for " + tokens[0]
              // +
              //                             "\n");
              // exit(0);
            }
            // vdwParams.insert({tokens[0], VdwParameters(std::stod(tokens[2]),
            // std::stod(tokens[3]))});
            // vdw14Params.insert(
            //     {tokens[0],
            //      VdwParameters(std::stod(tokens[5]), std::stod(tokens[6]))});
          } else {
            // std::cerr << "Extra tokens in line " << line << "\n";
            throw std::invalid_argument("Extra tokens in line " + line + "\n");
            exit(0);
          }
        }

      } // state : NONBONDED

      if (state == State::CMAP) {
        // std::cout << line << std::endl;
      }

      if (state == State::NBFIX) {
        // std::cout << "NBFIX : " << line << "\t";
        tokens = apo::split(line);
        if (tokens.size() >= 4) {
          if (tokens[0] > tokens[1])
            std::swap(tokens[0], tokens[1]);
          double emin = std::abs(std::stod(tokens[2]));
          double rmin = std::stod(tokens[3]);
          double emin14 = emin;
          double rmin14 = rmin;
          if (tokens.size() >= 5) {
            emin14 = std::stod(tokens[4]);
            rmin14 = std::stod(tokens[5]);
          }

          // std::cout << emin << " " << rmin << " " << emin14 << " " << rmin14
          //           << std::endl;
          NBFixParameters nbf{tokens[0], tokens[1], emin, rmin, emin14, rmin14};
          std::tuple<std::string, std::string> key{tokens[0], tokens[1]};
          // BondKey key{tokens[0], tokens[1]};
          //  nbfixParams[key] = nbf;
          nbfixParams.insert({key, nbf});
        }
      }
    }
  }
}
*/

void CharmmParameters::parseBondRecord(const std::vector<std::string> &tokens,
                                       const std::string &line,
                                       const std::string &fileName,
                                       const std::size_t lineNumber) {
  APOCHARMM_REQUIRE(tokens.size() == 4, ApoCharmmErrorCode::Runtime,
                    "Invalid BONDS parameter record in file \"" + fileName +
                        "\" at line " + std::to_string(lineNumber) + ": " +
                        line);

  std::string atom1 = tokens[0];
  std::string atom2 = tokens[1];
  if (atom1 > atom2)
    std::swap(atom1, atom2);

  const double kb = ParsePrmValue<double>(tokens[2], "kb", "BONDS", line,
                                          fileName, lineNumber);
  const double b0 = ParsePrmValue<double>(tokens[3], "b0", "BONDS", line,
                                          fileName, lineNumber);

  m_BondParams.insert({BondKey(atom1, atom2), BondValues(kb, b0)});

  return;
}

void CharmmParameters::parseAngleRecord(const std::vector<std::string> &tokens,
                                        const std::string &line,
                                        const std::string &fileName,
                                        const std::size_t lineNumber) {
  APOCHARMM_REQUIRE(
      (tokens.size() == 5) || (tokens.size() == 7), ApoCharmmErrorCode::Runtime,
      "Invalid ANGLES parameter record in file \"" + fileName + "\" at line " +
          std::to_string(lineNumber) + ": " + line);

  std::string atom1 = tokens[0];
  const std::string &atom2 = tokens[1];
  std::string atom3 = tokens[2];
  if (atom1 > atom3)
    std::swap(atom1, atom3);

  const double degreesToRadians = std::acos(-1.0) / 180.0;
  const double kTheta = ParsePrmValue<double>(tokens[3], "kTheta", "ANGLES",
                                              line, fileName, lineNumber);
  const double theta0 =
      degreesToRadians * ParsePrmValue<double>(tokens[4], "theta0", "ANGLES",
                                               line, fileName, lineNumber);

  double kub = 0.0;
  double s0 = 0.0;
  if (tokens.size() == 7) {
    kub = ParsePrmValue<double>(tokens[5], "kub", "ANGLES", line, fileName,
                                lineNumber);
    s0 = ParsePrmValue<double>(tokens[6], "s0", "ANGLES", line, fileName,
                               lineNumber);
  }

  const AngleKey key(atom1, atom2, atom3);
  m_AngleParams.insert({key, AngleValues(kTheta, theta0)});
  m_UreybParams.insert({key, BondValues(kub, s0)});

  return;
}

void CharmmParameters::parseDihedralRecord(
    const std::vector<std::string> &tokens, const std::string &line,
    const std::string &fileName, const std::size_t lineNumber) {
  APOCHARMM_REQUIRE(tokens.size() == 7, ApoCharmmErrorCode::Runtime,
                    "Invalid DIHEDRALS parameter record in file \"" + fileName +
                        "\" at line " + std::to_string(lineNumber) + ": " +
                        line);

  std::string atom1 = tokens[0];
  std::string atom2 = tokens[1];
  std::string atom3 = tokens[2];
  std::string atom4 = tokens[3];

  if (atom1 > atom4) {
    std::swap(atom1, atom4);
    std::swap(atom2, atom3);
  } else if ((atom1 == atom4) && (atom2 > atom3))
    std::swap(atom2, atom3);

  const double kChi = ParsePrmValue<double>(tokens[4], "kChi", "DIHEDRALS",
                                            line, fileName, lineNumber);
  const int multiplicity = ParsePrmValue<int>(
      tokens[5], "multiplicity", "DIHEDRALS", line, fileName, lineNumber);
  const double delta = ParsePrmValue<double>(tokens[6], "delta", "DIHEDRALS",
                                             line, fileName, lineNumber);

  const DihedralKey key(atom1, atom2, atom3, atom4);
  m_DihedralParams[key].push_back(DihedralValues(kChi, multiplicity, delta));

  return;
}

void CharmmParameters::parseImproperRecord(
    const std::vector<std::string> &tokens, const std::string &line,
    const std::string &fileName, const std::size_t lineNumber) {
  APOCHARMM_REQUIRE(tokens.size() == 7, ApoCharmmErrorCode::Runtime,
                    "Invalid IMPROPER parameter record in file \"" + fileName +
                        "\" at line " + std::to_string(lineNumber) + ": " +
                        line);

  std::string atom1 = tokens[0];
  std::string atom2 = tokens[1];
  std::string atom3 = tokens[2];
  std::string atom4 = tokens[3];

  if (atom1 > atom4) {
    std::swap(atom1, atom4);
    std::swap(atom2, atom3);
  } else if ((atom1 == atom4) && (atom2 > atom3))
    std::swap(atom2, atom3);

  const double kPsi = ParsePrmValue<double>(tokens[4], "kPsi", "IMPROPER", line,
                                            fileName, lineNumber);
  static_cast<void>(ParsePrmValue<double>(tokens[5], "ignored multiplicity",
                                          "IMPROPER", line, fileName,
                                          lineNumber));
  const double psi0 = ParsePrmValue<double>(tokens[6], "psi0", "IMPROPER", line,
                                            fileName, lineNumber);

  m_ImproperParams.insert(
      {DihedralKey(atom1, atom2, atom3, atom4), ImDihedralValues(kPsi, psi0)});

  return;
}

void CharmmParameters::parseNonbondedRecord(
    const std::vector<std::string> &tokens, const std::string &line,
    const std::string &fileName, const std::size_t lineNumber) {
  APOCHARMM_REQUIRE(
      (tokens.size() == 4) || (tokens.size() == 7), ApoCharmmErrorCode::Runtime,
      "Invalid NONBONDED parameter record in file \"" + fileName +
          "\" at line " + std::to_string(lineNumber) + ": " + line);

  static_cast<void>(ParsePrmValue<double>(
      tokens[1], "ignored value", "NONBONDED", line, fileName, lineNumber));
  const double epsilon = ParsePrmValue<double>(
      tokens[2], "epsilon", "NONBONDED", line, fileName, lineNumber);
  const double rmin_2 = ParsePrmValue<double>(tokens[3], "rmin/2", "NONBONDED",
                                              line, fileName, lineNumber);

  m_VdwParams.insert_or_assign(tokens[0], VdwParameters(epsilon, rmin_2));

  if (tokens.size() == 7) {
    static_cast<void>(ParsePrmValue<double>(
        tokens[4], "ignored 1-4", "NONBONDED", line, fileName, lineNumber));
    const double epsilon14 = ParsePrmValue<double>(
        tokens[5], "1-4 epsilon", "NONBONDED", line, fileName, lineNumber);
    const double rmin_2_14 = ParsePrmValue<double>(
        tokens[6], "1-4 rmin/2", "NONBONDED", line, fileName, lineNumber);

    m_Vdw14Params.insert_or_assign(tokens[0],
                                   VdwParameters(epsilon14, rmin_2_14));
  }

  return;
}

void CharmmParameters::parseNbfixRecord(const std::vector<std::string> &tokens,
                                        const std::string &line,
                                        const std::string &fileName,
                                        const std::size_t lineNumber) {
  APOCHARMM_REQUIRE(
      (tokens.size() == 4) || (tokens.size() == 6), ApoCharmmErrorCode::Runtime,
      "Invalid NBFIX parameter record in file \"" + fileName + "\" at line " +
          std::to_string(lineNumber) + ": " + line);

  std::string atom1 = tokens[0];
  std::string atom2 = tokens[1];
  if (atom1 > atom2)
    std::swap(atom1, atom2);

  const double emin = std::abs(ParsePrmValue<double>(
      tokens[2], "emin", "NBFIX", line, fileName, lineNumber));
  const double rmin = ParsePrmValue<double>(tokens[3], "rmin", "NBFIX", line,
                                            fileName, lineNumber);

  double emin14 = emin;
  double rmin14 = rmin;
  if (tokens.size() == 6) {
    emin14 = std::abs(ParsePrmValue<double>(tokens[4], "1-4 emin", "NBFIX",
                                            line, fileName, lineNumber));
    rmin14 = ParsePrmValue<double>(tokens[5], "1-4 rmin", "NBFIX", line,
                                   fileName, lineNumber);
  }

  const NBFixParameters parameters{atom1, atom2, emin, rmin, emin14, rmin14};
  const std::tuple<std::string, std::string> key{atom1, atom2};
  m_NbfixParams.insert({key, parameters});

  return;
}
