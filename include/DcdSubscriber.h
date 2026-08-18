// BEGINLICENSE
//
// This file is part of apoCHARMM, which is distributed under the BSD 3-clause
// license, as described in the LICENSE file in the top level directory of this
// project.
//
// Author: Samarjeet Prasad, James E. Gonzales II
//
// ENDLICENSE

#pragma once

#include "Subscriber.h"

#include <string>

/**
 * @brief Writes periodic coordinate snapshots in CHARMM DCD format.
 *
 * Construction immediately creates or truncates a binary output file. The first
 * successful @ref update writes a native-endian DCD header, title block, atom
 * count, and one frame. Later updates append frames and patch the header's
 * frame and step counters.
 *
 * Each frame contains one six-`double` unit-cell record followed by separate
 * X, Y, and Z blocks of single-precision coordinates in atom order. Coordinates
 * and box lengths use angstroms. Charges are not written. Coordinate output
 * transfers the complete single-precision coordinate/charge container from the
 * CUDA device to host memory and synchronizes the current CUDA device.
 *
 * @warning The binary representation uses native integer size, floating-point
 * representation, and byte order. The implementation performs no endian
 * conversion.
 * @warning Reporting frequency is cached by the integrator at subscription but
 * read again by DCD metadata generation. Changing it after subscription can
 * make callback spacing and recorded metadata disagree.
 * @warning A failed frame write is not transactional and can leave a partial or
 * internally inconsistent DCD file.
 *
 * @see Subscriber
 * @see subscriber
 */
class DcdSubscriber : public Subscriber {
public:
  /**
   * @brief Constructs a DCD writer with the default reporting interval.
   *
   * @param[in] fileName Output path copied by the subscriber. The value must be
   * nonempty; when it contains `/`, its nonempty parent path must exist.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `fileName` is empty or its checked
   * parent path does not exist.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if the
   * binary output file cannot be opened for writing.
   *
   * @post The file has been created or truncated, no DCD header has been
   * written, and `getReportFrequency() == 1000`.
   */
  DcdSubscriber(const std::string &fileName);

  /**
   * @brief Constructs a DCD writer with an explicit reporting interval.
   *
   * @param[in] fileName Output path copied by the subscriber. The value must be
   * nonempty; when it contains `/`, its nonempty parent path must exist.
   * @param[in] reportFrequency Positive, dimensionless number of propagated
   * steps between scheduled DCD frames.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if `reportFrequency` is not positive,
   * `fileName` is empty, or its checked parent path does not exist.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if the
   * binary output file cannot be opened for writing.
   *
   * @post The file has been created or truncated, no DCD header has been
   * written, and `getReportFrequency() == reportFrequency`.
   */
  DcdSubscriber(const std::string &fileName, const int reportFrequency);

  /**
   * @brief Destroys the writer and closes its output stream.
   *
   * Destruction is non-throwing.
   */
  ~DcdSubscriber(void) noexcept override = default;

public:
  /**
   * @brief Appends one DCD coordinate frame and updates header metadata.
   *
   * On the first update, the method writes the header before the unit-cell and
   * coordinate records. It then patches `NFILE` at byte offset 8 and `NSTEP` at
   * byte offset 20, flushes the stream, and increments the in-memory frame
   * count only after the final stream-state check succeeds.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotInitialized` if the logical file name is empty, no
   * context or integrator is attached, the stream is closed, the context does
   * not contain exactly three positive box lengths, or the atom count is not
   * positive.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if the next
   * frame count or derived step count exceeds signed `int`, or if header/frame
   * output or flushing fails.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if transferring
   * coordinates from device to host or synchronizing the current CUDA device
   * fails.
   *
   * @pre The context's single-precision coordinate/charge container has one
   * element per atom and valid device storage.
   * @post On success, one complete frame has been appended and the DCD header
   * reflects the new in-memory frame count.
   * @warning On failure, already-written bytes remain observable and the file
   * may require deletion rather than reuse.
   */
  void update(void) override;

  /**
   * @brief Opens the stored path as a truncating binary DCD output stream.
   *
   * An already-open stream is closed first. Successful opening resets the frame
   * count to zero and marks the header as unwritten.
   *
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::NotInitialized` if the stored file name is empty.
   * @throws ApoCharmmError With code
   * `ApoCharmmErrorCode::InvalidArgument` if the checked parent path does not
   * exist.
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if the file
   * cannot be opened for binary output.
   *
   * @post On success, prior file contents are truncated and the next update
   * writes a fresh header.
   * @warning A failed reopen can leave the previous stream closed while the
   * stored logical file name remains unchanged.
   */
  void openFile(void) override;

private:
  /**
   * @brief Writes the fixed DCD header, title, and atom-count records.
   *
   * The control record stores the current report frequency, context degree-of-
   * freedom count, integrator time step in picoseconds, an enabled unit-cell
   * flag, and CHARMM version value 49. It flushes and validates the stream
   * before marking the header written.
   *
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Runtime` if writing
   * or flushing the header fails.
   *
   * @pre The file is open and valid context and integrator attachments exist.
   */
  void writeHeader(void);

  /**
   * @brief Writes one six-value periodic-cell record.
   *
   * Elements 0, 2, and 5 receive the context's X, Y, and Z box lengths in
   * angstroms; the remaining `double` elements are zero.
   *
   * @pre The context contains exactly three valid box lengths and the stream is
   * open.
   */
  void writeXtalData(void);

  /**
   * @brief Transfers and writes the X, Y, and Z coordinate blocks.
   *
   * The complete context-owned `CudaContainer<float4>` is copied from device to
   * host and synchronized. Three Fortran-style records then store `float`
   * components in X-all-atoms, Y-all-atoms, Z-all-atoms order.
   *
   * @throws ApoCharmmError With code `ApoCharmmErrorCode::Cuda` if the transfer
   * or device-wide synchronization fails.
   *
   * @pre The context atom count is positive and agrees with both active
   * coordinate/charge mirror lengths.
   */
  void writeCoordData(void);

private:
  /**
   * @brief Stores the number of frames committed by successful updates.
   *
   * The value is encoded into the DCD `NFILE` field and participates in the
   * derived `NSTEP` value.
   */
  int m_NumFramesWritten;

  /** @brief Records whether the current open file already contains a header. */
  bool m_IsHeaderWritten;
};
