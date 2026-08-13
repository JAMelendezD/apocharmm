// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#include "Coordinates.h"

#include "ApoCharmmError.h"

#include <cstddef>
#include <string>
#include <vector_functions.h>

Coordinates::Coordinates(void)
    : m_NumAtoms(-1), m_CoordinatesDP(), m_CoordinatesSP() {}

Coordinates::Coordinates(const std::vector<double3> &coords) : Coordinates() {
  m_NumAtoms = static_cast<int>(coords.size());
  m_CoordinatesDP.resize(m_NumAtoms);
  m_CoordinatesSP.resize(m_NumAtoms);

  for (int i = 0; i < m_NumAtoms; i++) {
    m_CoordinatesDP[i] = make_double3(coords[i].x, coords[i].y, coords[i].z);
    m_CoordinatesSP[i] = make_float3(static_cast<float>(coords[i].x),
                                     static_cast<float>(coords[i].y),
                                     static_cast<float>(coords[i].z));
  }
}

Coordinates::Coordinates(const std::vector<float3> &coords) : Coordinates() {
  m_NumAtoms = static_cast<int>(coords.size());
  m_CoordinatesDP.resize(m_NumAtoms);
  m_CoordinatesSP.resize(m_NumAtoms);

  for (int i = 0; i < m_NumAtoms; i++) {
    m_CoordinatesDP[i] = make_double3(static_cast<double>(coords[i].x),
                                      static_cast<double>(coords[i].y),
                                      static_cast<double>(coords[i].z));
    m_CoordinatesSP[i] = make_float3(coords[i].x, coords[i].y, coords[i].z);
  }
}

Coordinates::Coordinates(const std::vector<std::vector<double>> &coords)
    : Coordinates() {
  m_NumAtoms = static_cast<int>(coords.size());
  m_CoordinatesDP.resize(m_NumAtoms);
  m_CoordinatesSP.resize(m_NumAtoms);

  for (int i = 0; i < m_NumAtoms; i++) {
    APOCHARMM_REQUIRE(coords[i].size() == 3,
                      ApoCharmmErrorCode::InvalidArgument,
                      "Coordinate entry " + std::to_string(i) +
                          " must contain exactly 3 values; observed " +
                          std::to_string(coords[i].size()));

    m_CoordinatesDP[i] = make_double3(coords[i][0], coords[i][1], coords[i][2]);
    m_CoordinatesSP[i] = make_float3(static_cast<float>(coords[i][0]),
                                     static_cast<float>(coords[i][1]),
                                     static_cast<float>(coords[i][2]));
  }
}

Coordinates::Coordinates(const std::vector<std::vector<float>> &coords)
    : Coordinates() {
  m_NumAtoms = static_cast<int>(coords.size());
  m_CoordinatesDP.resize(m_NumAtoms);
  m_CoordinatesSP.resize(m_NumAtoms);

  for (int i = 0; i < m_NumAtoms; i++) {
    APOCHARMM_REQUIRE(coords[i].size() == 3,
                      ApoCharmmErrorCode::InvalidArgument,
                      "Coordinate entry " + std::to_string(i) +
                          " must contain exactly 3 values; observed " +
                          std::to_string(coords[i].size()));

    m_CoordinatesDP[i] = make_double3(static_cast<double>(coords[i][0]),
                                      static_cast<double>(coords[i][1]),
                                      static_cast<double>(coords[i][2]));
    m_CoordinatesSP[i] = make_float3(coords[i][0], coords[i][1], coords[i][2]);
  }
}

void Coordinates::setNumAtoms(const int numAtoms) {
  APOCHARMM_REQUIRE(numAtoms >= 0, ApoCharmmErrorCode::InvalidArgument,
                    "Atom count must be non-negative; observed " +
                        std::to_string(numAtoms));

  m_NumAtoms = numAtoms;
  m_CoordinatesDP.resize(numAtoms);
  m_CoordinatesSP.resize(numAtoms);
  return;
}

const std::vector<double3> &Coordinates::getCoordinatesDP(void) const {
  return m_CoordinatesDP;
}

std::vector<double3> &Coordinates::getCoordinatesDP(void) {
  return m_CoordinatesDP;
}

const std::vector<float3> &Coordinates::getCoordinatesSP(void) const {
  return m_CoordinatesSP;
}

std::vector<float3> &Coordinates::getCoordinatesSP(void) {
  return m_CoordinatesSP;
}

int Coordinates::getNumAtoms(void) const { return m_NumAtoms; }
