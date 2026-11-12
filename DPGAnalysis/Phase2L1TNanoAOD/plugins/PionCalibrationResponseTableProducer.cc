// Calibration-only Nano response table. The matching, cone helpers and 10-bit
// rounding follow FastPUPPI ResponseNTuplizer.cc at d8521ad175c4a75157b7b30236e546f7d4fe5dfc.
// The selected prompt particle IDs are configured in Python. This supports the
// charged-pion barrel workflow and the mixed-particle HGCal workflow with the
// same matching definition.
// It deliberately does not implement the reference analyzer's tau or random-cone modes.
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/L1TParticleFlow/interface/PFCluster.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <map>
#include <vector>

namespace {
  struct SimpleObject {
    float pt, eta, phi, emf;
    SimpleObject(float apt, float aneta, float aphi, float aemf) : pt(apt), eta(aneta), phi(aphi), emf(aemf) {}
    bool operator<(const SimpleObject& other) const { return eta < other.eta; }
    bool operator<(const float& other) const { return eta < other; }
  };
  class MultiCollection {
  public:
    MultiCollection(const edm::ParameterSet& iConfig, const std::string& name, edm::ConsumesCollector&& coll)
        : name_(name), prop_(false), sel_("") {
      if (name.find("Ecal") != std::string::npos)
        prop_ = true;
      else if (name.find("Hcal") != std::string::npos)
        prop_ = true;
      else if (name.find("Calo") != std::string::npos)
        prop_ = true;
      const std::vector<edm::InputTag>& tags = iConfig.getParameter<std::vector<edm::InputTag>>(name);
      for (const auto& tag : tags)
        tokens_.push_back(coll.consumes<reco::CandidateView>(tag));
      if (iConfig.existsAs<std::string>(name + "_sel")) {
        sel_ = StringCutObjectSelector<reco::Candidate>(iConfig.getParameter<std::string>(name + "_sel"), true);
      }
    }
    const std::string& name() const { return name_; }
    bool prop() const { return prop_; }
    void get(const edm::Event& iEvent) {
      edm::Handle<reco::CandidateView> handle;
      for (const auto& token : tokens_) {
        iEvent.getByToken(token, handle);
        for (const reco::Candidate& c : *handle) {
          if (!sel_(c)) {
            continue;
          }
          const l1t::PFCluster* pfcluster = dynamic_cast<const l1t::PFCluster*>(&c);
          float emf = -1.;
          if (pfcluster) {
            auto encoding = pfcluster->encoding();
            if (encoding == l1t::PFCluster::HWEncoding::Had || encoding == l1t::PFCluster::HWEncoding::Em) {
              auto obj = pfcluster->caloDigiObj();
              if (auto digi = std::get_if<l1ct::HadCaloObj>(&obj)) {
                emf = digi->floatEmPt() / digi->floatPt();
              }
            }
          }
          objects_.emplace_back(c.pt(), c.eta(), c.phi(), emf);
        }
      }
      std::sort(objects_.begin(), objects_.end());
    }
    const std::vector<SimpleObject>& objects() const { return objects_; }
    void clear() { objects_.clear(); }

  private:
    std::string name_;
    bool prop_;
    std::vector<edm::EDGetTokenT<reco::CandidateView>> tokens_;
    StringCutObjectSelector<reco::Candidate> sel_;
    std::vector<SimpleObject> objects_;
  };
  class InCone {
  public:
    InCone(const std::vector<SimpleObject>& objects, float eta, float phi, float dr) {
      auto first =
          std::lower_bound(objects.begin(), objects.end(), eta - dr - 0.01f);  // small offset to avoid dealing with ==
      auto end = std::lower_bound(objects.begin(), objects.end(), eta + dr + 0.01f);
      float dr2 = dr * dr;
      sum04 = 0;
      for (auto it = first; it < end; ++it) {
        float mydr2 = ::deltaR2(eta, phi, it->eta, it->phi);
        if (mydr2 < dr2) {
          ptdr2.emplace_back(it->pt, mydr2);
          emfdr2.emplace_back(it->emf, mydr2);
        }
        if (mydr2 < 0.16f)
          sum04 += it->pt;
      }
    }
    float sum(float dr = 0.4) const {
      if (dr == 0.4f)
        return sum04;
      float dr2 = dr * dr;
      float mysum = 0;
      for (const auto& p : ptdr2) {
        if (p.second < dr2)
          mysum += p.first;
      }
      return mysum;
    }
    float emf_nearest() const {
      std::pair<float, float> best(0, 9999);
      for (const auto& emf : emfdr2) {
        if (emf.second < best.second)
          best = emf;
      }
      return best.first;
    }

  private:
    std::vector<std::pair<float, float>> ptdr2;
    std::vector<std::pair<float, float>> emfdr2;
    float sum04;
  };

}  // namespace
class PionCalibrationResponseTableProducer : public edm::stream::EDProducer<> {
public:
  explicit PionCalibrationResponseTableProducer(const edm::ParameterSet& cfg)
      : genJets_(consumes<std::vector<reco::GenJet>>(cfg.getParameter<edm::InputTag>("genJets"))),
        genParticles_(consumes<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"))),
        name_(cfg.getParameter<std::string>("name")),
        particlePdgIds_(cfg.getParameter<std::vector<int>>("particlePdgIds")) {
    const auto objects = cfg.getParameter<edm::ParameterSet>("objects");
    for (const auto& name : objects.getParameterNamesForType<std::vector<edm::InputTag>>())
      collections_.emplace_back(objects, name, consumesCollector());
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto& jets = event.get(genJets_);
    const auto& particles = event.get(genParticles_);
    std::vector<const reco::GenParticle*> prompts;
    for (const auto& particle : particles)
      if (particle.statusFlags().isPrompt())
        prompts.push_back(&particle);
    for (auto& collection : collections_)
      collection.get(event);

    // All columns share one row per selected gen-jet/particle match. Empty events
    // retain the same schema with zero rows. No fitting selection is imposed here.
    std::map<std::string, std::vector<float>> columns;
    for (const auto* field : {"mc_pt", "mc_eta", "mc_phi", "mc_caloeta", "mc_calophi", "mc_iso04"})
      columns[field] = {};
    for (const auto& collection : collections_)
      for (const auto* field : {"pt02", "emfbest"})
        columns[collection.name() + "_" + field] = {};
    std::vector<int> pdgIds;
    for (const auto& jet : jets) {
      const reco::GenParticle* match = nullptr;
      bool ambiguous = false;
      for (const auto* particle : prompts) {
        if (deltaR2(*particle, jet) < 0.16f) {
          if (match) {
            ambiguous = true;
            break;
          }
          match = particle;
        }
      }
      if (ambiguous || !match ||
          std::find(particlePdgIds_.begin(), particlePdgIds_.end(), std::abs(match->pdgId())) == particlePdgIds_.end())
        continue;
      const float pt = zip(match->pt());
      const float eta = zip(match->eta()), phi = zip(match->phi());
      const math::XYZTLorentzVector vertex(match->vx(), match->vy(), match->vz(), 0.);
      const auto propagated = l1tpf::propagateToCalo(match->p4(), vertex, match->charge(), 3.8112);
      const float caloeta = zip(propagated.first), calophi = zip(propagated.second);
      pdgIds.push_back(std::abs(match->pdgId()));
      columns["mc_pt"].push_back(pt);
      columns["mc_eta"].push_back(eta);
      columns["mc_phi"].push_back(phi);
      columns["mc_caloeta"].push_back(caloeta);
      columns["mc_calophi"].push_back(calophi);
      columns["mc_iso04"].push_back(zip(jet.pt() / pt - 1));
      for (const auto& collection : collections_) {
        const InCone cone(
            collection.objects(), collection.prop() ? caloeta : eta, collection.prop() ? calophi : phi, 0.8);
        columns[collection.name() + "_pt02"].push_back(zip(cone.sum(0.2)));
        columns[collection.name() + "_emfbest"].push_back(zip(cone.emf_nearest()));
      }
    }
    auto table = std::make_unique<nanoaod::FlatTable>(pdgIds.size(), name_, false, false);
    table->addColumn<int>("mc_id", pdgIds, "Absolute generator PDG ID; prompt particles matched uniquely to gen jets");
    for (const auto& column : columns)
      table->addColumn<float>(column.first, column.second, columnDocumentation(column.first), -1);
    event.put(std::move(table));
    for (auto& collection : collections_)
      collection.clear();
  }

private:
  static std::string columnDocumentation(const std::string& name) {
    if (name == "mc_pt")
      return "Matched generator-particle transverse momentum [GeV]";
    if (name == "mc_eta")
      return "Matched generator-particle pseudorapidity";
    if (name == "mc_phi")
      return "Matched generator-particle azimuth [rad]";
    if (name == "mc_caloeta")
      return "Generator-particle pseudorapidity propagated to the calorimeter";
    if (name == "mc_calophi")
      return "Generator-particle azimuth propagated to the calorimeter [rad]";
    if (name == "mc_iso04")
      return "Generator jet pt divided by rounded pion pt, minus one";
    if (name.find("_pt02") != std::string::npos)
      return "Sum of cluster pt within deltaR < 0.2 of the generator direction [GeV]";
    return "EM pt fraction of the nearest cluster within deltaR < 0.8; zero if no cluster";
  }

  static float zip(float value) { return MiniFloatConverter::reduceMantissaToNbitsRounding<10>(value); }
  edm::EDGetTokenT<std::vector<reco::GenJet>> genJets_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_;
  std::string name_;
  std::vector<int> particlePdgIds_;
  std::vector<MultiCollection> collections_;
};
DEFINE_FWK_MODULE(PionCalibrationResponseTableProducer);
