// -*- C++ -*-
//
// Package: L1CaloTrigger
// Class: Phase2L1TJetSeedEmulatorProducer
//
/**\class Phase2L1TJetSeedEmulatorProducer Phase2L1TJetSeedEmulatorProducer.cc L1Trigger/L1CaloTrigger/plugin/Phase2L1TJetSeedEmulatorProducer.cc
*/

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/L1TParticleFlow/interface/puppi.h"
#include "DataFormats/L1TParticleFlow/interface/gt_datatypes.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/common/bitonic_hybrid_sort_ref.h"

#include <cmath>

#include <algorithm>
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1TJetSeedEmulator.h"

class Phase2L1TJetSeedEmulatorProducer : public edm::one::EDProducer<> {
public:
  explicit Phase2L1TJetSeedEmulatorProducer(const edm::ParameterSet&);
  ~Phase2L1TJetSeedEmulatorProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  void convertEDMToHW(const l1t::PFCandidateCollection&, std::vector<l1ct::PuppiObj>&);

  edm::EDGetTokenT<l1t::PFCandidateCollection> inputCollectionTag_;
  
  bool debug;
  size_t nBinsEta;
  unsigned int nBinsPhi;
  unsigned int jetIEtaSize;
  unsigned int jetIPhiSize;
  bool trimmedGrid;
  double seedPtThreshold;
  // Eta and phi edges of input PF regions
  std::vector<double> etaRegionEdges;
  std::vector<double> phiRegionEdges;
  // Maximum number of candidates per input PF region
  unsigned int maxInputsPerRegion;
  Phase2L1TJetSeedEmulator emulator;
  std::string outputCollectionName;

};

Phase2L1TJetSeedEmulatorProducer::Phase2L1TJetSeedEmulatorProducer(const edm::ParameterSet& iConfig)
  : inputCollectionTag_{
      consumes<l1t::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("inputCollectionTag"))},
      debug(iConfig.getParameter<bool>("debug")),
      nBinsEta(iConfig.getParameter<unsigned int>("nBinsEta")),
      nBinsPhi(iConfig.getParameter<unsigned int>("nBinsPhi")),
      jetIEtaSize(iConfig.getParameter<unsigned int>("jetIEtaSize")),
      jetIPhiSize(iConfig.getParameter<unsigned int>("jetIPhiSize")),
      trimmedGrid(iConfig.getParameter<bool>("trimmedGrid")),
      seedPtThreshold(iConfig.getParameter<double>("seedPtThreshold")),
      etaRegionEdges(iConfig.getParameter<std::vector<double>>("etaRegions")),
      phiRegionEdges(iConfig.getParameter<std::vector<double>>("phiRegions")),
      maxInputsPerRegion(iConfig.getParameter<unsigned int>("maxInputsPerRegion")),
      emulator(debug, nBinsEta, nBinsPhi, jetIEtaSize, jetIPhiSize, trimmedGrid, 
          seedPtThreshold, etaRegionEdges, phiRegionEdges, maxInputsPerRegion),
      outputCollectionName(iConfig.getParameter<std::string>("outputCollectionName")) {

  produces<l1t::PFCandidateCollection>(outputCollectionName);
}

Phase2L1TJetSeedEmulatorProducer::~Phase2L1TJetSeedEmulatorProducer() {}


void Phase2L1TJetSeedEmulatorProducer::convertEDMToHW(const l1t::PFCandidateCollection& inputCollection, std::vector<l1ct::PuppiObj>& puppiObjects) {
  puppiObjects.reserve(inputCollection.size());
  for (const auto& candidate : inputCollection) {
    l1ct::PuppiObj puppiObj;
    puppiObj.initFromBits(candidate.encodedPuppi64());
    puppiObjects.push_back(puppiObj);
  }
}

void Phase2L1TJetSeedEmulatorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<l1t::PFCandidateCollection> inputCollectionHandle;
  iEvent.getByToken(inputCollectionTag_, inputCollectionHandle);

  std::vector<l1ct::PuppiObj> puppiObjects;
  convertEDMToHW(*inputCollectionHandle, puppiObjects);

  l1t::PFCandidateCollection sortedSeeds = emulator.emulateEvent(puppiObjects);

  std::unique_ptr<l1t::PFCandidateCollection> outputSeedsCollection(new l1t::PFCandidateCollection);
  outputSeedsCollection->swap(sortedSeeds);

  iEvent.put(std::move(outputSeedsCollection), outputCollectionName);

  return;
}



void Phase2L1TJetSeedEmulatorProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<bool>("debug", false);
  desc.add<edm::InputTag>("inputCollectionTag", edm::InputTag("l1tLayer1", "Puppi"));
  desc.add<unsigned int>("nBinsEta", 72);
  desc.add<unsigned int>("nBinsPhi", 72);
  desc.add<unsigned int>("jetIEtaSize", 9);
  desc.add<unsigned int>("jetIPhiSize", 9);
  desc.add<bool>("trimmedGrid", true);
  desc.add<double>("seedPtThreshold", 1);
  desc.add<std::string>("outputCollectionName", "histoJetSeeds9x9trimmed");
  desc.add<std::vector<double>>("etaRegions", { -3, -2.5, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2.5, 3 });
  desc.add<std::vector<double>>("phiRegions", { -3.15, -2.45, -1.75, -1.05, -0.35, 0.35, 1.05, 1.75, 2.45, 3.15 });
  desc.add<unsigned int>("maxInputsPerRegion", 18);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(Phase2L1TJetSeedEmulatorProducer);
