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
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
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
  

  edm::EDGetTokenT<edm::View<reco::Candidate>> inputCollectionTag_;
  
  bool debug;
  std::vector<double> etaBinning;
  size_t nBinsEta;
  unsigned int nBinsPhi;
  unsigned int jetIEtaSize;
  unsigned int jetIPhiSize;
  bool trimmedGrid;
  double seedPtThreshold;
  double ptlsb;
  double philsb;
  double etalsb;
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
          consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("inputCollectionTag"))},
      debug(iConfig.getParameter<bool>("debug")),
      etaBinning(iConfig.getParameter<std::vector<double>>("etaBinning")),
      nBinsEta(etaBinning.size() - 1),
      nBinsPhi(iConfig.getParameter<unsigned int>("nBinsPhi")),
      jetIEtaSize(iConfig.getParameter<unsigned int>("jetIEtaSize")),
      jetIPhiSize(iConfig.getParameter<unsigned int>("jetIPhiSize")),
      trimmedGrid(iConfig.getParameter<bool>("trimmedGrid")),
      seedPtThreshold(iConfig.getParameter<double>("seedPtThreshold")),
      ptlsb(iConfig.getParameter<double>("ptlsb")),
      philsb(iConfig.getParameter<double>("philsb")),
      etalsb(iConfig.getParameter<double>("etalsb")),
      etaRegionEdges(iConfig.getParameter<std::vector<double>>("etaRegions")),
      phiRegionEdges(iConfig.getParameter<std::vector<double>>("phiRegions")),
      maxInputsPerRegion(iConfig.getParameter<unsigned int>("maxInputsPerRegion")),
      emulator(debug, etaBinning, nBinsPhi, jetIEtaSize, jetIPhiSize, trimmedGrid, 
               seedPtThreshold, ptlsb, philsb, etalsb, etaRegionEdges, phiRegionEdges, maxInputsPerRegion),
      outputCollectionName(iConfig.getParameter<std::string>("outputCollectionName")) {

    produces<l1t::PFCandidateCollection>(outputCollectionName);
}

Phase2L1TJetSeedEmulatorProducer::~Phase2L1TJetSeedEmulatorProducer() {}


void Phase2L1TJetSeedEmulatorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Candidate>> inputCollectionHandle;
  iEvent.getByToken(inputCollectionTag_, inputCollectionHandle);

  std::unique_ptr<l1t::PFCandidateCollection> jetSeedsCollection(new l1t::PFCandidateCollection);

  l1t::PFCandidateCollection sortedSeeds = emulator.emulateEvent( inputCollectionHandle );

  jetSeedsCollection->swap(sortedSeeds);

  iEvent.put(std::move(jetSeedsCollection), outputCollectionName );

  return;
}



void Phase2L1TJetSeedEmulatorProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<bool>("debug", false);
  desc.add<edm::InputTag>("inputCollectionTag", edm::InputTag("l1tLayer1", "Puppi"));
  desc.add<std::vector<double>>("etaBinning", { -3 , -2.91491519955 , -2.83201206065 , -2.74910892175 , -2.66620578285 , -2.58330264395 , -2.5 , -2.41491519955 , -2.33201206065 , -2.24910892175 , -2.16620578285 , -2.08330264395 , -2.00039950505 , -1.91749636615 , -1.83459322725 , -1.75169008835 , -1.66878694945 , -1.58588381055 , -1.5 , -1.41491519955 , -1.33201206065 , -1.24910892175 , -1.16620578285 , -1.08330264395 , -1.0 , -0.91491519955 , -0.83201206065 , -0.74910892175 , -0.66620578285 , -0.58330264395 , -0.5 , -0.41491519955 , -0.33201206065 , -0.24910892175 , -0.16620578285 , -0.08330264395 , 0 , 0.08508480045 , 0.16798793935 , 0.25089107825 , 0.33379421715 , 0.41669735605 , 0.5 , 0.58508480045 , 0.66798793935 , 0.75089107825 , 0.83379421715 , 0.91669735605 , 1 , 1.08508480045 , 1.16798793935 , 1.25089107825 , 1.33379421715 , 1.41669735605 , 1.5 , 1.58508480045 , 1.66798793935 , 1.75089107825 , 1.83379421715 , 1.91669735605 , 1.99960049495 , 2.08250363385 , 2.16540677275 , 2.24830991165 , 2.33121305055 , 2.41411618945 , 2.5 , 2.58508480045 , 2.66798793935 , 2.75089107825 , 2.83379421715 , 2.91669735605 , 3});
  desc.add<unsigned int>("nBinsPhi", 72);
  desc.add<unsigned int>("jetIEtaSize", 9);
  desc.add<unsigned int>("jetIPhiSize", 9);
  desc.add<bool>("trimmedGrid", true);
  desc.add<double>("seedPtThreshold", 1);
  desc.add<double>("ptlsb", 0.25),
  desc.add<double>("philsb", 0.0043633231),
  desc.add<double>("etalsb", 0.0043633231),
  desc.add<string>("outputCollectionName", "histoJetSeeds9x9trimmed");
  desc.add<std::vector<double>>("etaRegions", { -3, -2.5, -1.5, -1.0, -0.5, 0, 0.5, 1, 1.5, 2.5, 3 });
  desc.add<std::vector<double>>("phiRegions", { -3.15, -2.45, -1.75, -1.05, -0.35, 0.35, 1.05, 1.75, 2.45, 3.15 });
  desc.add<unsigned int>("maxInputsPerRegion", 18);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(Phase2L1TJetSeedEmulatorProducer);
