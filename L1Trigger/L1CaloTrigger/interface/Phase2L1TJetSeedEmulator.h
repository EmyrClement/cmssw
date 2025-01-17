#ifndef L1Trigger_L1CaloTrigger_Phase2L1TJetSeedEmulator_h
#define L1Trigger_L1CaloTrigger_Phase2L1TJetSeedEmulator_h
// -*- C++ -*-
//
// Package:     L1Trigger/L1CaloTrigger
// Class  :     Phase2L1TJetSeedEmulator
//
/**\class Phase2L1TJetSeedEmulator Phase2L1TJetSeedEmulator.h "Phase2L1TJetSeedEmulator.h"

 Description: HistoSeededCone Seed Finding

*/
//
// Original Author:  Dharmender
// Created:  Tue, 03 Dec 2024 15:29:22 GMT
//

// system include files
#include "DataFormats/L1TParticleFlow/interface/puppi.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/common/bitonic_hybrid_sort_ref.h"

#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>

class Phase2L1TJetSeedEmulator {
public:
  Phase2L1TJetSeedEmulator(bool debug, unsigned int nBinsEta, unsigned int nBinsPhi, unsigned int jetIEtaSize, unsigned int jetIPhiSize, bool trimmedGrid, double seedPtThreshold, std::vector<double> etaRegionEdges, std::vector<double> phiRegionEdges ,unsigned int maxInputsPerRegion );

  std::vector<l1ct::PuppiObj> emulateEvent(const std::vector<l1ct::PuppiObj>& puppiObjects);

  std::vector<l1ct::PuppiObj> findSeeds(float seedThreshold) const;
  float getBinContent(int iEta, int iPhi) const;
  bool trimBin(int etaIndex, int phiIndex) const;
  void sortSeeds(const std::vector<l1ct::PuppiObj>& unsortedSeeds, std::vector<l1ct::PuppiObj>& sortedSeeds);

  std::pair<double, double> regionEtaPhiLowEdges(unsigned int regionIndex) const;
  std::pair<double, double> regionEtaPhiUpEdges(unsigned int regionIndex) const;
  std::pair<unsigned, unsigned> regionEtaPhiBinOffset(unsigned int regionIndex) const;
  std::pair<unsigned, unsigned> getCandidateBin(const l1ct::glbeta_t glbEta, const l1ct::glbphi_t glbPhi, const unsigned int regionIndex) const;

  template <typename T>
  void swap(T& a, T& b);

  template <typename T>
  void compAndSwap(std::vector<T>& a, unsigned int i, unsigned int j, bool dir=false);

  template <typename T>
  void hybridBitonicMergeRef(std::vector<T>& a, int N, int low, bool dir);

  template <typename T>
  void hybridBitonicSortRef(std::vector<T>& a, int N, int low, bool dir);

  template <typename T>
  void hybrid_bitonic_sort_and_crop_ref(unsigned int nIn, unsigned int nOut, const std::vector<T>& in, std::vector<T>& out);

  void fillHistogram(std::vector<std::vector<l1ct::pt_t>>& histogram, const std::vector<l1ct::PuppiObj>& puppis, unsigned int regionIndex);

  unsigned int getRegionIndex(unsigned int phiRegion, unsigned int etaRegion) const;

  std::vector<std::vector<l1ct::PuppiObj>> prepareInputsIntoRegions(const std::vector<l1ct::PuppiObj>& puppiObjects);

private:
  bool debug_;

  unsigned int nBinsEta_;
  unsigned int nBinsPhi_;
  unsigned int jetIEtaSize_;
  unsigned int jetIPhiSize_;
  bool trimmedGrid_;
  double seedPtThreshold_;
  std::vector<double> etaRegionEdges_;
  std::vector<double> phiRegionEdges_;
  unsigned int maxInputsPerRegion_;
  double etaBinLSB_;
  double phiBinLSB_;
  int etaBinSize_;
  int phiBinSize_;
  unsigned int nBinsPhiRegion_; // New data member
  std::vector<std::vector<l1ct::pt_t>> histogram_;

  // Constants used by seed sort 
  static constexpr unsigned int nEtaRegions_ = 4;
  static constexpr unsigned int nInputsPerSortModule_ = 18;
  static constexpr unsigned int nOutputSeedsPerEtaRegion_ = 4;
  static constexpr unsigned int nOutputSeedsToGT_ = 12;
};

template <typename T>
void Phase2L1TJetSeedEmulator::swap(T& a, T& b) {
  T temp = a;
  a = b;
  b = temp;
}

template <typename T>
void Phase2L1TJetSeedEmulator::compAndSwap(std::vector<T>& a, unsigned int i, unsigned int j, bool dir) {
  if (i >= a.size() || j >= a.size() || i == j) return;

  if (dir) {
    if (a[j] < a[i]) std::swap(a[i], a[j]);
  } else {
    if (a[i] < a[j]) std::swap(a[i], a[j]);
  }
}

template <typename T>
void Phase2L1TJetSeedEmulator::hybridBitonicMergeRef(std::vector<T>& a, int N, int low, bool dir) {
  int k = hybridBitonicSortUtils::PowerOf2LessThan(N);
  int k2 = N - k;

  if (N > 1) {
    for (int i = low; i < low + k; i++) {
      if (i + k < low + N) compAndSwap(a, i, i + k, dir);
    }
    if (N > 2) {
      hybridBitonicMergeRef(a, k, low, dir);
      hybridBitonicMergeRef(a, k2, low + k, dir);
    }
  }
}

template <typename T>
void Phase2L1TJetSeedEmulator::hybridBitonicSortRef(std::vector<T>& a, int N, int low, bool dir) {
  if (N > 1) {
    int lowerSize = N / 2;
    int upperSize = N - lowerSize;
    hybridBitonicSortRef(a, lowerSize, low, !dir);
    hybridBitonicSortRef(a, upperSize, low + lowerSize, dir);
    hybridBitonicMergeRef(a, N, low, dir);
  }
}

template <typename T>
void Phase2L1TJetSeedEmulator::hybrid_bitonic_sort_and_crop_ref(unsigned int nIn, unsigned int nOut, const std::vector<T>& in, std::vector<T>& out) {
  std::vector<T> work = in;
  hybridBitonicSortRef(work, nIn, 0, false);

  out.resize(nOut);
  for (unsigned int i = 0; i < nOut; ++i) {
    out[i] = work[i];
  }
}

#endif