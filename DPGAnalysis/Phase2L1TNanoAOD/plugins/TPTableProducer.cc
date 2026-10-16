#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using SimpleEcalEBTPFlatTableProducer = SimpleFlatTableProducer<EcalEBTriggerPrimitiveDigi>;
DEFINE_FWK_MODULE(SimpleEcalEBTPFlatTableProducer);

#if __has_include("DataFormats/EcalDigi/interface/EcalEBPhase2TriggerPrimitiveDigi.h")
using SimpleEcalEBPhase2TPFlatTableProducer = SimpleFlatTableProducer<EcalEBPhase2TriggerPrimitiveDigi>;
DEFINE_FWK_MODULE(SimpleEcalEBPhase2TPFlatTableProducer);
#endif

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/L1THGCal/interface/HGCalTower.h"

using SimpleHcalTPFlatTableProducer = SimpleFlatTableProducer<HcalTriggerPrimitiveDigi>;
using SimpleTPTrackFlatTableProducer = SimpleFlatTableProducer<TTTrack<Ref_Phase2TrackerDigi_>>;
using SimpleTPHGCalClusterFlatTableProducer = BXVectorSimpleFlatTableProducer<l1t::HGCalMulticluster>;
using SimpleTPHGCalTowerFlatTableProducer = BXVectorSimpleFlatTableProducer<l1t::HGCalTower>;
DEFINE_FWK_MODULE(SimpleHcalTPFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTPTrackFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTPHGCalClusterFlatTableProducer);
DEFINE_FWK_MODULE(SimpleTPHGCalTowerFlatTableProducer);
