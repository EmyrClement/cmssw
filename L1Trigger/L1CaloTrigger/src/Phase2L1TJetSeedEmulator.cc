// -*- C++ -*-
//
// Package:     L1Trigger/L1CaloTrigger
// Class  :     Phase2L1TJetSeedEmulator
//
// Implementation:
//     [Notes on implementation]
//
// Original Author:  Dharmender
//         Created:  Tue, 03 Dec 2024 15:29:22 GMT
//

// system include files

// user include files
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1TJetSeedEmulator.h"
//
// constructors and destructor
//
Phase2L1TJetSeedEmulator::Phase2L1TJetSeedEmulator(bool debug, unsigned int nBinsEta, unsigned int nBinsPhi, unsigned int jetIEtaSize, unsigned int jetIPhiSize, bool trimmedGrid, double seedPtThreshold, std::vector<double> etaRegionEdges, std::vector<double> phiRegionEdges ,unsigned int maxInputsPerRegion) 
  : debug_(debug),
    nBinsEta_(nBinsEta),
    nBinsPhi_(nBinsPhi),
    jetIEtaSize_(jetIEtaSize),
    jetIPhiSize_(jetIPhiSize),
    trimmedGrid_(trimmedGrid),
    seedPtThreshold_(seedPtThreshold),
    etaRegionEdges_(etaRegionEdges),
    phiRegionEdges_(phiRegionEdges),
    maxInputsPerRegion_(maxInputsPerRegion),
    etaBinLSB_((etaRegionEdges.back() - etaRegionEdges.front()) / nBinsEta),
    phiBinLSB_((phiRegionEdges.back() - phiRegionEdges.front()) / nBinsPhi),
    etaBinSize_(l1ct::Scales::makeGlbEta(etaBinLSB_)),
    phiBinSize_(l1ct::Scales::makeGlbPhi(phiBinLSB_)),
    nBinsPhiRegion_( nBinsPhi_ / ( phiRegionEdges_.size() - 1 ) ),
    histogram_(nBinsEta, std::vector<l1ct::pt_t>(nBinsPhi, 0)) { 
}

bool Phase2L1TJetSeedEmulator::trimBin(const int etaIndex, const int phiIndex) const {
  int etaHalfSize = jetIEtaSize_ / 2;
  int phiHalfSize = jetIPhiSize_ / 2;

  if (etaIndex == -etaHalfSize || etaIndex == etaHalfSize) {
    if (phiIndex <= -phiHalfSize + 1 || phiIndex >= phiHalfSize - 1) {
      return true;
    }
  } else if (etaIndex == -etaHalfSize + 1 || etaIndex == etaHalfSize - 1) {
    if (phiIndex == -phiHalfSize || phiIndex == phiHalfSize) {
      return true;
    }
  }

  return false;
}

// Phase2L1TJetSeedEmulator::~Phase2L1TJetSeedEmulator() {}

//
// member functions
//

std::vector<l1ct::PuppiObj> Phase2L1TJetSeedEmulator::emulateEvent(const std::vector<std::vector<l1ct::PuppiObj>>& puppiObjects2D, const std::vector<std::pair<double, double>>& regionLowEdges) {
  // Resetting histogram
  for (auto& row : histogram_) {
    std::fill(row.begin(), row.end(), 0);
  }
  // histogramming the data
  for (unsigned int iInputRegion = 0; iInputRegion < puppiObjects2D.size(); ++iInputRegion) {
    if (puppiObjects2D[iInputRegion].empty()) {
      continue;
    }
    double etaLowEdge = regionLowEdges[iInputRegion].first;
    double phiLowEdge = regionLowEdges[iInputRegion].second;
    fillHistogram(histogram_, puppiObjects2D[iInputRegion], etaLowEdge, phiLowEdge);
  }

  // find the seeds
  const auto& seedsVector = findSeeds(seedPtThreshold_);

  // sort by pt
  std::vector<l1ct::PuppiObj> sortedSeeds;
  sortSeeds(seedsVector, sortedSeeds);

  return sortedSeeds;
}

float Phase2L1TJetSeedEmulator::getBinContent(int iEta, int iPhi) const {
  int nBinsEta = histogram_.size();
  int nBinsPhi = histogram_[0].size();
  while (iPhi < 0) {
    iPhi += nBinsPhi;
  }
  while (iPhi >= nBinsPhi) {
    iPhi -= nBinsPhi;
  }
  if (iEta < 0 || iEta >= nBinsEta) {
    return 0;
  }
  return histogram_[iEta][iPhi];
}

std::vector<l1ct::PuppiObj> Phase2L1TJetSeedEmulator::findSeeds(float seedThreshold) const {
  int nBinsX = histogram_.size();
  int nBinsY = histogram_[0].size();

  std::vector<l1ct::PuppiObj> seeds;

  int etaHalfSize = (int)jetIEtaSize_ / 2;
  int phiHalfSize = (int)jetIPhiSize_ / 2;

  for (int iPhi = 0; iPhi < nBinsY; iPhi++) {
    for (int iEta = 0; iEta < nBinsX; iEta++) {
      float centralPt = histogram_[iEta][iPhi];
      if (centralPt < seedThreshold)
        continue;

      bool isLocalMaximum = true;
      for (int etaIndex = -etaHalfSize; etaIndex <= etaHalfSize; etaIndex++) {
        for (int phiIndex = -phiHalfSize; phiIndex <= phiHalfSize; phiIndex++) {
          if (trimmedGrid_) {
            if (trimBin(etaIndex, phiIndex))
              continue;
          }

          if ((etaIndex == 0) && (phiIndex == 0))
            continue;
          if (etaIndex > 0) {
            isLocalMaximum = ((isLocalMaximum) && (centralPt > getBinContent(iEta + etaIndex, iPhi + phiIndex)));
          } else if ( etaIndex < 0 ) {
            isLocalMaximum = ((isLocalMaximum) && (centralPt >= getBinContent(iEta + etaIndex, iPhi + phiIndex)));
          }
          else {
            if ( phiIndex > 0 ) {
              isLocalMaximum = ((isLocalMaximum) && (centralPt > getBinContent(iEta + etaIndex, iPhi + phiIndex)));
            }
            else {
              isLocalMaximum = ((isLocalMaximum) && (centralPt >= getBinContent(iEta + etaIndex, iPhi + phiIndex)));
            }
          }
        }
      }

      if (isLocalMaximum) {
        l1ct::PuppiObj seed;
        seed.hwPt = centralPt;
        seed.hwEta = l1ct::Scales::makeGlbEta(etaRegionEdges_.front() + (iEta + 0.5) * etaBinLSB_);
        seed.hwPhi = l1ct::Scales::makeGlbPhi(phiRegionEdges_.front() + (iPhi + 0.5) * phiBinLSB_);
        seeds.emplace_back(seed);
      }
    }
  }
  return seeds;
}

void Phase2L1TJetSeedEmulator::sortSeeds(const std::vector<l1ct::PuppiObj>& unsortedSeeds, std::vector<l1ct::PuppiObj>& sortedSeeds) {
  // Get seeds into the regions and time ordering seen in firmware
  std::vector<std::vector<std::vector<std::vector<l1ct::PuppiObj>>>> seedsPerEtaPhiRegions(
    nEtaRegions_, std::vector<std::vector<std::vector<l1ct::PuppiObj>>>(
      2, std::vector<std::vector<l1ct::PuppiObj>>(
        nInputsPerSortModule_, std::vector<l1ct::PuppiObj>())));
        for (const auto& seed : unsortedSeeds) {
          unsigned int etaRegion = (seed.hwEta + l1ct::Scales::makeGlbEta(3)) / l1ct::Scales::makeGlbEta(1.5);
          unsigned int seedPhiBin = (seed.hwPhi + l1ct::Scales::makeGlbPhi(M_PI)) / phiBinSize_;
          unsigned int phiRegion = ((seedPhiBin) % 4) / 2;

          if (etaRegion >= nEtaRegions_ || phiRegion >= 2 || seedPhiBin / 4 >= nInputsPerSortModule_) {
            continue;
          }

          seedsPerEtaPhiRegions[etaRegion][phiRegion][seedPhiBin / 4].push_back(seed);
        }

        // Rotate to first phi region found in firmware
        for (unsigned iEtaRegion = 0; iEtaRegion < nEtaRegions_; ++iEtaRegion) {
          for (unsigned iPhiRegion = 0; iPhiRegion < 2; ++iPhiRegion) {
            std::rotate(seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin(), seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin() + 8, seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].end());
          }
        }

  // Push seeds in first phi bin to back, as these are found last after receiving all bins (i.e. handling of phi wrap-around)
  for (unsigned iEtaRegion = 0; iEtaRegion < nEtaRegions_; ++iEtaRegion) {
    for (unsigned iPhiRegion = 0; iPhiRegion < 2; ++iPhiRegion) {
      std::rotate(seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin(), seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin() + 1, seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].end());
    }
  }

  std::vector<l1ct::PuppiObj> sortedSeedsAllEta;
  for (unsigned iEtaRegion = 0; iEtaRegion < nEtaRegions_; ++iEtaRegion) {
    std::vector<l1ct::PuppiObj> sortedSeedsInEtaRegion;
    for (unsigned iPhiRegion = 0; iPhiRegion < 2; ++iPhiRegion) {
      std::vector<l1ct::PuppiObj> sortedSeeds(4, l1ct::PuppiObj());
      for (unsigned int iInputClock = 0; iInputClock < nInputsPerSortModule_; ++iInputClock) {
        // Sort input seeds
        std::vector<l1ct::PuppiObj> inputSeeds = seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion][iInputClock];
        // First by eta
        std::sort(inputSeeds.begin(), inputSeeds.end(), [](l1ct::PuppiObj seed1, l1ct::PuppiObj seed2) {
          return seed1.hwEta < seed2.hwEta;
        });
        inputSeeds.resize(nOutputSeedsPerEtaRegion_);
        hybrid_bitonic_sort_and_crop_ref(4, 4, inputSeeds, inputSeeds);

        // Add to list of top 4 seeds so far
        // Merge with top 4 seeds so far, and sort
        sortedSeeds.insert(sortedSeeds.end(), inputSeeds.begin(), inputSeeds.end());
        std::reverse(sortedSeeds.begin(), sortedSeeds.begin() + nOutputSeedsPerEtaRegion_);
        for (int i = 0; i < 4; i++) {
          compAndSwap(sortedSeeds, i, i + 4, 0);
        }

        sortedSeeds.resize(nOutputSeedsPerEtaRegion_);
        std::reverse(sortedSeeds.begin(), sortedSeeds.end());
        compAndSwap(sortedSeeds, 0, 2);
        compAndSwap(sortedSeeds, 1, 3);
        //---
        compAndSwap(sortedSeeds, 0, 1);
        compAndSwap(sortedSeeds, 2, 3);
      }

      if (iPhiRegion % 2 == 0) {
        sortedSeedsInEtaRegion.insert(sortedSeedsInEtaRegion.end(), sortedSeeds.rbegin(), sortedSeeds.rend());
      } else {
        sortedSeedsInEtaRegion.insert(sortedSeedsInEtaRegion.end(), sortedSeeds.begin(), sortedSeeds.end());
      }
    }
    // Sort 8 seeds in each eta region
    std::reverse(sortedSeedsInEtaRegion.begin(), sortedSeedsInEtaRegion.end());
    hybridBitonicMergeRef(sortedSeedsInEtaRegion, nOutputSeedsPerEtaRegion_ * 2, 0, false);

    if (iEtaRegion % 2 == 0) {
      sortedSeedsAllEta.insert(sortedSeedsAllEta.end(), sortedSeedsInEtaRegion.rbegin(), sortedSeedsInEtaRegion.rend());
    } else {
      sortedSeedsAllEta.insert(sortedSeedsAllEta.end(), sortedSeedsInEtaRegion.begin(), sortedSeedsInEtaRegion.end());
    }
  }
  hybridBitonicMergeRef(sortedSeedsAllEta, nOutputSeedsPerEtaRegion_ * 2 * 2, 0, false);
  hybridBitonicMergeRef(sortedSeedsAllEta, nOutputSeedsPerEtaRegion_ * 2 * 2, nOutputSeedsPerEtaRegion_ * 2 * 2, false);
  std::reverse(sortedSeedsAllEta.begin(), sortedSeedsAllEta.begin() + nOutputSeedsPerEtaRegion_ * 2 * 2);

  for (unsigned int iJet = 0; iJet < nOutputSeedsPerEtaRegion_ * 2 * 2 - nOutputSeedsToGT_; ++iJet) {
    sortedSeedsAllEta.erase(sortedSeedsAllEta.begin());
    sortedSeedsAllEta.erase(sortedSeedsAllEta.end() - 1);
  }

  hybridBitonicMergeRef(sortedSeedsAllEta, nOutputSeedsToGT_ * 2, 0, false);
  sortedSeedsAllEta.resize(nOutputSeedsToGT_);
  unsigned int nSeedsGT0 = 0;
  for (const auto& iJet : sortedSeedsAllEta) {
    if (iJet.hwPt > 0) {
      sortedSeeds.push_back(iJet);
      ++nSeedsGT0;
    }
  }
}

std::pair<unsigned, unsigned> Phase2L1TJetSeedEmulator::regionEtaPhiBinOffset(double etaLowEdge, double phiLowEdge) const {
  float phiRegionWidth = abs(phiRegionEdges_.at(0) - phiRegionEdges_.at(1));
  float phiBinOffset = (phiLowEdge - phiRegionEdges_.front()) / phiRegionWidth * nBinsPhiRegion_;

  float etaBinOffset = (etaLowEdge + 3) / 0.5 * 6;

  return std::pair<unsigned, unsigned>{phiBinOffset, etaBinOffset};
}

std::pair<unsigned, unsigned> Phase2L1TJetSeedEmulator::getCandidateBin(const l1ct::glbeta_t glbEta, const l1ct::glbphi_t glbPhi, double etaLowEdge, double phiLowEdge) const {
  l1ct::glbeta_t etaOffset = l1ct::Scales::makeGlbEta(etaLowEdge);
  l1ct::glbphi_t phiOffset = l1ct::Scales::makeGlbPhi(phiLowEdge);

  // Debug printout
  std::cout << "getCandidateBin calculations:" << std::endl;
  std::cout << "  glbEta: " << glbEta << ", glbPhi: " << glbPhi << std::endl;
  std::cout << "  etaOffset: " << etaOffset << ", phiOffset: " << phiOffset << std::endl;

  int etaBin = (glbEta - etaOffset) / etaBinSize_ + 1;
  int phiBin = (glbPhi - phiOffset) / phiBinSize_ + 1;

  constexpr int nBinsEtaRegionWithTrack = 12;
  constexpr int nBinsEtaRegionEverywhereElse = 6;
  if (etaLowEdge == -2.5 || etaLowEdge == 1.5) {
    if (etaBin >= nBinsEtaRegionWithTrack) etaBin = nBinsEtaRegionWithTrack;
  } else if (etaBin >= nBinsEtaRegionEverywhereElse) {
    etaBin = nBinsEtaRegionEverywhereElse;
  }
  if (phiBin >= int(nBinsPhiRegion_)) phiBin = nBinsPhiRegion_;

  // Hopefully temporary fix for handling candidates with phi=pi
  if (glbPhi == 720) phiBin = 1;

  std::pair<unsigned, unsigned> binOffsets = regionEtaPhiBinOffset(etaLowEdge, phiLowEdge);

  // Debug printout
  std::cout << "  etaBin: " << etaBin << ", phiBin: " << phiBin << std::endl;
  std::cout << "  binOffsets.first: " << binOffsets.first << ", binOffsets.second: " << binOffsets.second << std::endl;

  return std::pair<unsigned, unsigned>{phiBin + binOffsets.first - 1, etaBin + binOffsets.second - 1};
}

void Phase2L1TJetSeedEmulator::fillHistogram(std::vector<std::vector<l1ct::pt_t>>& histogram, const std::vector<l1ct::PuppiObj>& puppis, double etaLowEdge, double phiLowEdge) {
  for (const auto& puppi : puppis) {
    std::cout << "Binning candidate : " << puppi.hwPt << " " << puppi.hwEta << " " << puppi.hwPhi << std::endl;
    auto binEtaPhi = getCandidateBin(puppi.hwEta, puppi.hwPhi, etaLowEdge, phiLowEdge);
    std::cout << "Got bins : " << binEtaPhi.first << " " << binEtaPhi.second << std::endl;
    histogram[binEtaPhi.second][binEtaPhi.first] += puppi.hwPt;
  }
}

unsigned int Phase2L1TJetSeedEmulator::getRegionIndex(const unsigned int phiRegion, const unsigned int etaRegion) const {
  return etaRegion * (phiRegionEdges_.size() - 1) + phiRegion;
}

