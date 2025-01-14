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
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1TParticleFlow/interface/puppi.h"
#include "DataFormats/L1TParticleFlow/interface/gt_datatypes.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/common/bitonic_hybrid_sort_ref.h"

#include "TH2F.h"

#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>

class Phase2L1TJetSeedEmulator {
public:
  Phase2L1TJetSeedEmulator(bool debug, std::vector<double> etaBinning, unsigned int nBinsPhi, unsigned int jetIEtaSize, unsigned int jetIPhiSize, bool trimmedGrid, double seedPtThreshold, double ptlsb, double philsb, double etalsb, std::vector<double> etaRegionEdges, std::vector<double> phiRegionEdges ,unsigned int maxInputsPerRegion);

  template <class Handle>
  l1t::PFCandidateCollection emulateEvent( Handle triggerPrimitives );

  l1t::PFCandidateCollection findSeeds(float seedThreshold) const;
  float getTowerEnergy(int iEta, int iPhi) const;
  bool trimTower(int etaIndex, int phiIndex) const;
  void sortSeeds(const l1t::PFCandidateCollection unsortedJets, l1t::PFCandidateCollection& sortedJets);

  std::pair<double, double> regionEtaPhiLowEdges(unsigned int regionIndex) const;
  std::pair<double, double> regionEtaPhiUpEdges(unsigned int regionIndex) const;
  std::pair<unsigned, unsigned> regionEtaPhiBinOffset(unsigned int regionIndex) const;
  std::pair<unsigned, unsigned> getCandidateBin(float eta, float phi, unsigned int regionIndex) const;

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

  template <class Container>
  void fillCaloGrid(TH2F& caloGrid, const Container& triggerPrimitives, unsigned int regionIndex);

  unsigned int getRegionIndex(unsigned int phiRegion, unsigned int etaRegion) const;

  template <class Handle>
  std::vector<std::vector<edm::Ptr<reco::Candidate>>> prepareInputsIntoRegions(const Handle& triggerPrimitives);

private:
  bool debug_;
  std::unique_ptr<TH2F> caloGrid_;

  std::vector<double> etaBinning_;
  size_t nBinsEta_;
  unsigned int nBinsPhi_;
  unsigned int jetIEtaSize_;
  unsigned int jetIPhiSize_;
  bool trimmedGrid_;
  double seedPtThreshold_;
  double ptlsb_;
  double philsb_;
  double etalsb_;
  std::vector<double> etaRegionEdges_;
  std::vector<double> phiRegionEdges_;
  unsigned int maxInputsPerRegion_;
};

// Template implementations
template <class Handle>
l1t::PFCandidateCollection Phase2L1TJetSeedEmulator::emulateEvent( Handle triggerPrimitives  ) {
  // sort inputs into PF regions
  std::vector<std::vector<reco::CandidatePtr>> inputsInRegions = prepareInputsIntoRegions<Handle>(triggerPrimitives);

  // histogramming the data
  caloGrid_->Reset();
  for (unsigned int iInputRegion = 0; iInputRegion < inputsInRegions.size(); ++iInputRegion) {
    fillCaloGrid<>(*(caloGrid_), inputsInRegions[iInputRegion], iInputRegion);
  }

  // int nBinsX = caloGrid_->GetNbinsX();
  // int nBinsY = caloGrid_->GetNbinsY();
  // for (int iPhi = 1; iPhi <= nBinsY; iPhi++)
  // {
  //   std::cout << "iPhi " << iPhi - 1 << " " << caloGrid_->GetYaxis()->GetBinCenter(iPhi) << " " << l1gt::phi_t(caloGrid_->GetYaxis()->GetBinCenter(iPhi) / l1gt::Scales::ETAPHI_LSB ) << ": ";
  //   for (int iEta = 1; iEta <= nBinsX; iEta++)
  //   {
  //     std::cout <<caloGrid_->GetBinContent(iEta, iPhi) << " ";
  //   }
  //   std::cout << std::endl;
  // }

  // find the seeds
  const auto& seedsVector = findSeeds(seedPtThreshold_);  // seedPtThreshold = 5

  // sort by pt
  l1t::PFCandidateCollection sortedSeeds;
  sortSeeds( seedsVector, sortedSeeds );
  return sortedSeeds;
}


template <typename T>
void Phase2L1TJetSeedEmulator::swap(T& a, T& b) {
  T temp = a;
  a = b;
  b = temp;
}

template <class Container>
void Phase2L1TJetSeedEmulator::fillCaloGrid(TH2F& caloGrid, const Container& triggerPrimitives, unsigned int regionIndex) {
  for (const auto& primitive : triggerPrimitives) {
    auto binEtaPhi = getCandidateBin(primitive->eta(), primitive->phi(), regionIndex);
    unsigned int globalBin = caloGrid.GetBin(binEtaPhi.second, binEtaPhi.first);
    caloGrid.AddBinContent(globalBin, float(l1ct::pt_t(primitive->pt())));
  }
}

template <typename T>
void Phase2L1TJetSeedEmulator::compAndSwap(std::vector<T>& a, unsigned int i, unsigned int j, bool dir) {
  if (i >= a.size() || j >= a.size() || i == j) return;

  if (dir) {
    if (a[j].pt() < a[i].pt()) std::swap(a[i], a[j]);
  } else {
    if (a[i].pt() < a[j].pt()) std::swap(a[i], a[j]);
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



template <class Handle>
std::vector<std::vector<edm::Ptr<reco::Candidate>>> Phase2L1TJetSeedEmulator::prepareInputsIntoRegions(const Handle& triggerPrimitives) {
  std::vector<std::vector<reco::CandidatePtr>> inputsInRegions{etaRegionEdges_.size() * (phiRegionEdges_.size() - 1)};

  for (unsigned int i = 0; i < triggerPrimitives->size(); ++i) {
    reco::CandidatePtr tp(triggerPrimitives, i);

    if (
      tp->phi() < phiRegionEdges_.front() || tp->phi() >= phiRegionEdges_.back() ||
        tp->eta() < etaRegionEdges_.front() || tp->eta() >= etaRegionEdges_.back())
      continue;

    // Which phi region does this tp belong to
    auto it_phi = phiRegionEdges_.begin();
    auto tp_phi = tp->phi();

    it_phi = std::upper_bound(phiRegionEdges_.begin(), phiRegionEdges_.end(), tp_phi) - 1;
    if ( l1ct::Scales::makeGlbPhi( *(it_phi+1) ) == l1ct::Scales::makeGlbPhi( tp_phi ) ) {
      it_phi += 1;
    }
    // Which eta region does this tp belong to
    auto it_eta = etaRegionEdges_.begin();
    it_eta = std::upper_bound(etaRegionEdges_.begin(), etaRegionEdges_.end(), tp->eta()) - 1;
    if ( l1ct::Scales::makeGlbEta( *(it_eta+1) ) == l1ct::Scales::makeGlbEta( tp->eta() ) ) {
      it_eta += 1;
    }


    if (it_phi != phiRegionEdges_.end() && it_eta != etaRegionEdges_.end()) {
      auto phiRegion = it_phi - phiRegionEdges_.begin();
      auto etaRegion = it_eta - etaRegionEdges_.begin();
      inputsInRegions[getRegionIndex(phiRegion, etaRegion)].emplace_back(tp);
    }
  }

  // Truncate number of inputs in each pf region
  for (auto& inputs : inputsInRegions) {
    if (inputs.size() > maxInputsPerRegion_) {
      inputs.resize(maxInputsPerRegion_);
    }
  }

  return inputsInRegions;
}

#endif