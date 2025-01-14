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

l1t::PFCandidateCollection Phase2L1TJetSeedEmulator::findSeeds(float seedThreshold) const {
  int nBinsX = histogram_.size();
  int nBinsY = histogram_[0].size();

  l1t::PFCandidateCollection seeds;

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
        l1t::PFCandidate p;
        reco::Candidate::PolarLorentzVector pfVector;

        const float etaLSB = 1.5 / 18;
        double etaBinCentre = -3 + (iEta+0.5)*etaLSB;

        const float phiLSB = 2. * M_PI / 72;
        double phiBinCentre = -M_PI + ( iPhi+0.5 ) * phiLSB;

        pfVector.SetPt(centralPt);
        pfVector.SetPhi(phiBinCentre);
        pfVector.SetEta(etaBinCentre);
        p.setP4( pfVector );

        l1ct::PuppiObj puppiObj;
        puppiObj.hwPt = l1ct::Scales::makePtFromFloat( centralPt );
        puppiObj.hwEta = l1ct::Scales::makeGlbEta( etaBinCentre );
        puppiObj.hwPhi = l1ct::Scales::makeGlbPhi( phiBinCentre );
        p.setEncodedPuppi64( puppiObj.pack().to_uint64() );

        seeds.emplace_back(p);
      }
    }
  }
  return seeds;
}

void Phase2L1TJetSeedEmulator::sortSeeds(const l1t::PFCandidateCollection unsortedSeeds, l1t::PFCandidateCollection& sortedSeeds ) {

  const unsigned int nEtaRegions = 4;
  const unsigned int nInputsPerSortModule = 18;
  const unsigned int nOutputSeedsPerEtaRegion = 4;
  const unsigned int nOutputSeedsToGT = 12;

  // unsigned int nUnsortedSeeds = unsortedSeeds.size();
  // Get seeds into the regions and time ordering seen in firmware
  std::vector< std::vector< std::vector< l1t::PFCandidateCollection > > > seedsPerEtaPhiRegions( 
    nEtaRegions, std::vector< std::vector< l1t::PFCandidateCollection > > (
      2, std::vector< l1t::PFCandidateCollection > (
        nInputsPerSortModule, l1t::PFCandidateCollection() ) ) );


  for ( const auto& seed : unsortedSeeds ) { 
    unsigned int etaRegion = (seed.eta()+3)/1.5;
    // unsigned int seedEtaBin = floor( ( seed.eta() + (2 - 1.0*etaRegion) * 1.5 ) / 0.0833 );
    unsigned int seedPhiBin = floor( ( seed.phi() + M_PI ) / 0.0875 );
    unsigned int phiRegion = ( ( seedPhiBin ) % 4 ) / 2;
    seedsPerEtaPhiRegions[etaRegion][phiRegion][seedPhiBin/4].push_back(seed);
  }

  // Rotate to first phi region found in firmware
  for ( unsigned iEtaRegion = 0; iEtaRegion < nEtaRegions; ++iEtaRegion ) {
    for ( unsigned iPhiRegion = 0; iPhiRegion < 2; ++ iPhiRegion ) {
      std::rotate( seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin(), seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin()+8, seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].end() );
    }
  }

  // Push seeds in first phi bin to back, as these are found last after receiving all bins (i.e. handling of phi wrap-around)
  for ( unsigned iEtaRegion = 0; iEtaRegion < nEtaRegions; ++iEtaRegion ) {
    for ( unsigned iPhiRegion = 0; iPhiRegion < 2; ++ iPhiRegion ) {
      std::rotate( seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin(), seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].begin()+1, seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion].end() );
    }
  }
  std::vector< l1t::PFCandidate > sortedSeedsAllEta;
  for ( unsigned iEtaRegion = 0; iEtaRegion < nEtaRegions; ++iEtaRegion ) {
    std::vector<l1t::PFCandidate > sortedSeedsInEtaRegion;
    for ( unsigned iPhiRegion = 0; iPhiRegion < 2; ++ iPhiRegion ) {
      std::vector<l1t::PFCandidate > sortedSeeds( 4, l1t::PFCandidate() );
      for ( unsigned int iInputClock = 0; iInputClock < nInputsPerSortModule; ++iInputClock ) {

        // Sort input seeds
        l1t::PFCandidateCollection inputSeeds = seedsPerEtaPhiRegions[iEtaRegion][iPhiRegion][iInputClock];
        // First by eta
        std::sort(inputSeeds.begin(), inputSeeds.end(), [](l1t::PFCandidate seed1, l1t::PFCandidate seed2) {
          return seed1.eta() < seed2.eta();
        });
        inputSeeds.resize(nOutputSeedsPerEtaRegion);
        hybrid_bitonic_sort_and_crop_ref(4,4,inputSeeds,inputSeeds);

        // Add to list of top 4 seeds so far
        // Merge with top 4 seeds so far, and sort
        sortedSeeds.insert( sortedSeeds.end(), inputSeeds.begin(), inputSeeds.end() );
        std::reverse(sortedSeeds.begin(),sortedSeeds.begin()+nOutputSeedsPerEtaRegion);
        for (int i = 0; i < 4; i++) {
            compAndSwap(sortedSeeds, i, i + 4, 0);
        }

        sortedSeeds.resize(nOutputSeedsPerEtaRegion);
        std::reverse(sortedSeeds.begin(),sortedSeeds.end());
        compAndSwap(sortedSeeds, 0, 2); 
        compAndSwap(sortedSeeds, 1, 3); 
        //---
        compAndSwap(sortedSeeds, 0, 1); 
        compAndSwap(sortedSeeds, 2, 3); 
      }

      if ( iPhiRegion % 2 == 0 ) {
        sortedSeedsInEtaRegion.insert( sortedSeedsInEtaRegion.end(), sortedSeeds.rbegin(), sortedSeeds.rend() );
      }
      else {
        sortedSeedsInEtaRegion.insert( sortedSeedsInEtaRegion.end(), sortedSeeds.begin(), sortedSeeds.end() );
      }
    }
    // Sort 8 seeds in each eta region
    // std::cout << "8 seeds in one of the regions, before merge" << std::endl;
    std::reverse(sortedSeedsInEtaRegion.begin(),sortedSeedsInEtaRegion.end());
    hybridBitonicMergeRef(sortedSeedsInEtaRegion,nOutputSeedsPerEtaRegion*2,0,false);

    if ( iEtaRegion % 2 == 0 ) {
      sortedSeedsAllEta.insert(sortedSeedsAllEta.end(), sortedSeedsInEtaRegion.rbegin(), sortedSeedsInEtaRegion.rend() );
    }
    else {
      sortedSeedsAllEta.insert(sortedSeedsAllEta.end(), sortedSeedsInEtaRegion.begin(), sortedSeedsInEtaRegion.end() );
    }
  }
  hybridBitonicMergeRef(sortedSeedsAllEta,nOutputSeedsPerEtaRegion*2*2,0,false);
  hybridBitonicMergeRef(sortedSeedsAllEta,nOutputSeedsPerEtaRegion*2*2,nOutputSeedsPerEtaRegion*2*2,false);
  std::reverse(sortedSeedsAllEta.begin(),sortedSeedsAllEta.begin()+nOutputSeedsPerEtaRegion*2*2);

  for ( unsigned int iJet = 0; iJet < nOutputSeedsPerEtaRegion*2*2 - nOutputSeedsToGT; ++iJet ) {
    sortedSeedsAllEta.erase(sortedSeedsAllEta.begin());
    sortedSeedsAllEta.erase(sortedSeedsAllEta.end()-1);
  }

  hybridBitonicMergeRef(sortedSeedsAllEta,nOutputSeedsToGT*2,0,false);
  sortedSeedsAllEta.resize(nOutputSeedsToGT);
  unsigned int nSeedsGT0=0;
  for ( const auto& iJet : sortedSeedsAllEta ) {
    if ( iJet.pt() > 0 ) {
      sortedSeeds.push_back( iJet );
      ++nSeedsGT0;
    }
  }
}

  // std::sort(seeds.begin(), seeds.end(), [](const l1t::PFCandidate& a, const l1t::PFCandidate& b) {    //sorting seeds by pt --should we use the regionised approach?
  //   return a.pt() > b.pt();
  // });
std::pair<double, double> Phase2L1TJetSeedEmulator::regionEtaPhiLowEdges(const unsigned int regionIndex) const {
  unsigned int phiRegion = regionIndex % (phiRegionEdges_.size() - 1);
  unsigned int etaRegion = (regionIndex - phiRegion) / (phiRegionEdges_.size() - 1);
  return std::pair<double, double>{phiRegionEdges_.at(phiRegion), etaRegionEdges_.at(etaRegion)};
}

std::pair<unsigned, unsigned> Phase2L1TJetSeedEmulator::regionEtaPhiBinOffset(const unsigned int regionIndex) const {
  unsigned int phiRegion = regionIndex % (phiRegionEdges_.size() - 1);
  unsigned int etaRegion = (regionIndex - phiRegion) / (phiRegionEdges_.size() - 1);

  float etaBinOffset = ( 3 + etaRegionEdges_.at(etaRegion) ) / 0.5 * 6;
  float phiRegionWidth = abs(phiRegionEdges_.at(0) - phiRegionEdges_.at(1) );

  float nBinsPhiRegion = round(phiRegionWidth/(2*pi/72));
  float phiBinOffset = ( -1.0 * phiRegionEdges_.front() + phiRegionEdges_.at(phiRegion) ) / phiRegionWidth * nBinsPhiRegion;
  return std::pair<unsigned, unsigned>{phiBinOffset, etaBinOffset};
}

std::pair<double, double> Phase2L1TJetSeedEmulator::regionEtaPhiUpEdges(const unsigned int regionIndex) const {
  unsigned int phiRegion = regionIndex % (phiRegionEdges_.size() - 1);
  unsigned int etaRegion = (regionIndex - phiRegion) / (phiRegionEdges_.size() - 1);
  if (phiRegion == phiRegionEdges_.size() - 1) {
    return std::pair<double, double>{phiRegionEdges_.at(phiRegion), etaRegionEdges_.at(etaRegion + 1)};
  } else if (etaRegion == etaRegionEdges_.size() - 1) {
    return std::pair<double, double>{phiRegionEdges_.at(phiRegion + 1), etaRegionEdges_.at(etaRegion)};
  }

  return std::pair<double, double>{phiRegionEdges_.at(phiRegion + 1), etaRegionEdges_.at(etaRegion + 1)};
}

std::pair<unsigned, unsigned> Phase2L1TJetSeedEmulator::getCandidateBin(const l1ct::glbeta_t glbEta, const l1ct::glbphi_t glbPhi, const unsigned int regionIndex) const {
  std::pair<double, double> regionLowEdges = regionEtaPhiLowEdges(regionIndex);
  l1ct::glbeta_t etaOffset = l1ct::Scales::makeGlbEta(regionLowEdges.second);
  l1ct::glbphi_t phiOffset = l1ct::Scales::makeGlbPhi(regionLowEdges.first);

  int etaBin = (glbEta - etaOffset) / 19 + 1;
  int phiBin = (glbPhi - phiOffset) / 20 + 1;

  if (regionLowEdges.second == -2.5 || regionLowEdges.second == 1.5) {
    if (etaBin >= 12) etaBin = 12;
  } else if (etaBin >= 6) {
    etaBin = 6;
  }
  if (phiBin >= 8) phiBin = 8;

  std::pair<unsigned, unsigned> binOffsets = regionEtaPhiBinOffset(regionIndex);

  return std::pair<unsigned, unsigned>{phiBin + binOffsets.first - 1, etaBin + binOffsets.second - 1};
}

unsigned int Phase2L1TJetSeedEmulator::getRegionIndex(const unsigned int phiRegion, const unsigned int etaRegion) const {
  return etaRegion * (phiRegionEdges_.size() - 1) + phiRegion;
}

l1t::PFCandidateCollection Phase2L1TJetSeedEmulator::emulateEvent(const std::vector<l1ct::PuppiObj>& puppiObjects) {
  // sort inputs into PF regions
  std::vector<std::vector<l1ct::PuppiObj>> inputsInRegions = prepareInputsIntoRegions(puppiObjects);

  // Resetting histogram
  for (auto& row : histogram_) {
    std::fill(row.begin(), row.end(), 0);
  }
  // histogramming the data
  for (unsigned int iInputRegion = 0; iInputRegion < inputsInRegions.size(); ++iInputRegion) {
    fillHistogram(histogram_, inputsInRegions[iInputRegion], iInputRegion);
  }

  // find the seeds
  const auto& seedsVector = findSeeds(seedPtThreshold_);

  // sort by pt
  l1t::PFCandidateCollection sortedSeeds;
  sortSeeds(seedsVector, sortedSeeds);
  return sortedSeeds;
}

void Phase2L1TJetSeedEmulator::fillHistogram(std::vector<std::vector<l1ct::pt_t>>& histogram, const std::vector<l1ct::PuppiObj>& puppis, unsigned int regionIndex) {
  for (const auto& puppi : puppis) {
    auto binEtaPhi = getCandidateBin(puppi.hwEta, puppi.hwPhi, regionIndex);
    histogram[binEtaPhi.second][binEtaPhi.first] += puppi.hwPt;
  }
}

std::vector<std::vector<l1ct::PuppiObj>> Phase2L1TJetSeedEmulator::prepareInputsIntoRegions(const std::vector<l1ct::PuppiObj>& puppiObjects) {
  std::vector<std::vector<l1ct::PuppiObj>> inputsInRegions(etaRegionEdges_.size() * (phiRegionEdges_.size() - 1));

  for (const auto& tp : puppiObjects) {
    if (tp.hwPhi < l1ct::Scales::makeGlbPhi(phiRegionEdges_.front()) || tp.hwPhi >= l1ct::Scales::makeGlbPhi(phiRegionEdges_.back()) ||
        tp.hwEta < l1ct::Scales::makeGlbEta(etaRegionEdges_.front()) || tp.hwEta >= l1ct::Scales::makeGlbEta(etaRegionEdges_.back())) {
      continue;
    }

    // Which phi region does this tp belong to
    auto it_phi = phiRegionEdges_.begin();
    it_phi = std::upper_bound(phiRegionEdges_.begin(), phiRegionEdges_.end(), l1ct::Scales::floatPhi(tp.hwPhi)) - 1;
    if (l1ct::Scales::makeGlbPhi(*(it_phi + 1)) == tp.hwPhi) {
      it_phi += 1;
    }

    // Which eta region does this tp belong to
    auto it_eta = etaRegionEdges_.begin();
    it_eta = std::upper_bound(etaRegionEdges_.begin(), etaRegionEdges_.end(), l1ct::Scales::floatEta(tp.hwEta)) - 1;
    if (l1ct::Scales::makeGlbEta(*(it_eta + 1)) == tp.hwEta) {
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

