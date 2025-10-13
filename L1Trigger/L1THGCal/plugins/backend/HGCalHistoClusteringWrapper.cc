#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "L1Trigger/L1THGCal/interface/HGCalAlgoWrapperBase.h"

#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_semiemulator/Stage2.hh"

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerModuleDetId.h"

#include "L1Trigger/L1THGCal/interface/backend/HGCalShowerShape.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"

class HGCalHistoClusteringWrapper : public HGCalHistoClusteringWrapperBase {
public:
  HGCalHistoClusteringWrapper(const edm::ParameterSet& conf);
  ~HGCalHistoClusteringWrapper() override {}

  void configure(
      const std::tuple<const HGCalTriggerGeometryBase* const, const edm::ParameterSet&, const unsigned int, const int>&
          configuration) override;

  void process(const std::vector<edm::Ptr<l1t::HGCalCluster>>& inputClusters,
               std::pair<l1t::HGCalMulticlusterBxCollection&, l1t::HGCalClusterBxCollection&>&
                   outputMulticlustersAndRejectedClusters) const override;

private:
  void convertCMSSWInputs(const std::vector<edm::Ptr<l1t::HGCalCluster>>& clustersPtrs,
                          std::vector<TPGTCBits>& clusters_SA) const;
  void convertAlgorithmOutputs(std::vector<TPGCluster>& clusters_SA_out,
                               l1t::HGCalMulticlusterBxCollection& multiClusters_out,
                               const std::vector<edm::Ptr<l1t::HGCalCluster>>& tcPtrs) const;

  void makeClusters(const std::vector<TPGTCBits>& triggerCells_in_SA,
                    std::vector<TPGCluster>& clusters_SA_out) const;

  void setGeometry(const HGCalTriggerGeometryBase* const geom) { 
    triggerTools_.setGeometry(geom);
    shape_.setGeometry(geom);
  }

  double rotatePhiToSectorZero(const double phi, const unsigned sector) const;
  double rotatePhiFromSectorZero(const double phi, const unsigned sector) const;

  HGCalTriggerTools triggerTools_;

  l1thgcfirmware::ClusterAlgoConfig theConfiguration_;

  mutable TPGStage2Emulation::Stage2 theAlgo_;

  edm::ESHandle<HGCalTriggerGeometryBase> triggerGeometry_;
  HGCalShowerShape shape_;

  edm::FileInPath muEtaLut_;
  edm::FileInPath sigmaEtaLut_;
};

HGCalHistoClusteringWrapper::HGCalHistoClusteringWrapper(const edm::ParameterSet& conf)
    : HGCalHistoClusteringWrapperBase(conf),
      theConfiguration_(),
      theAlgo_(conf.getParameterSet("layer2FwClusteringParameters").getParameter<double>("sideLength")* sqrt(3.0)),
      muEtaLut_("L1Trigger/L1THGCal/data/mean_eta_LUT.csv"),
      sigmaEtaLut_("L1Trigger/L1THGCal/data/sigma_eta_LUT.csv") {}

void HGCalHistoClusteringWrapper::convertCMSSWInputs(const std::vector<edm::Ptr<l1t::HGCalCluster>>& clustersPtrs,
                                                     std::vector<TPGTCBits>& clusters_SA) const {

  for (size_t idx = 0; idx < clustersPtrs.size(); ++idx) {
    const auto& clusterPtr = clustersPtrs[idx];
    TPGTCFloats tcFloat;
    tcFloat.setZero();
    tcFloat.setROverZPhiF(clusterPtr->position().x() / std::abs( clusterPtr->position().z() ),
                          clusterPtr->position().y() / std::abs( clusterPtr->position().z() ),
                          theConfiguration_.sector());
    tcFloat.setEnergyGeV(clusterPtr->pt());
    DetId id(clusterPtr->detId());
    // Semi emulator expects global layer, not trigger layer
    tcFloat.setLayer(triggerTools_.layerWithOffset(id));

    // unsigned layer = triggerTools_.layerWithOffset(id);
    // if ( layer <= 13 ) layer *= 2;
    // else layer += 13;

    // // unsigned layer = 0;
    // unsigned det = DetId(id).det();
    // // unsigned subdet = 0;

    // // if (det == DetId::HGCalTrigger) {
    // //   subdet = HGCalTriggerDetId(id).subdet();
    // //   if (subdet == HGCalTriggerSubdetector::HGCalEETrigger) {
    // //     layer = HGCalTriggerDetId(id).layer();
    // //   } else if (subdet == HGCalTriggerSubdetector::HGCalHSiTrigger) {
    // //     layer = heOffset_ + HGCalTriggerDetId(id).layer();
    // //   }
    // // } else if (det == DetId::HGCalHSc) {
    // //   layer = heOffset_ + HGCScintillatorDetId(id).layer();
    // // } else if (det == DetId::Forward) {
    // //   subdet = HGCalTriggerModuleDetId(id).triggerSubdetId();
    // //   if (subdet == HGCalTriggerSubdetector::HGCalEETrigger) {
    // //     layer = HGCalTriggerModuleDetId(id).layer();
    // //   } else if (subdet == HGCalTriggerSubdetector::HGCalHSiTrigger ||
    // //             subdet == HGCalTriggerSubdetector::HGCalHScTrigger) {
    // //     layer = HGCalDetId(id).layer();
    // //   }
    // // }
    // std::cout << "Trigger cell det : " << det
    //           << ", layer : " << layer << " " << triggerTools_.layerWithOffset(id) << std::endl;
    // tcFloat.setLayer(layer);
    tcFloat.setCMSSWIndex(idx);

    // if ( clusterPtr->position().z() > 0 ) continue;
    // if ( clusterPtr->phi() < -1.96 ||  clusterPtr->phi() > -1.26 ) continue;
    // if ( clusterPtr->eta() < -3. ||  clusterPtr->eta() > -2.03 ) continue;
    // if ( triggerTools_.layerWithOffset(id) <= 5 ) continue;  // Only use layers 27-40
    // if ( clusterPtr->pt() > 5 ) {
    // std::cout << "Input TC: "
    //           << "Energy: " << clusterPtr->pt() << ", "
    //           << "Eta: " << clusterPtr->eta() << ", "
    //           << "Phi: " << clusterPtr->phi() << ", "
    //           << "Layer: " << triggerTools_.layerWithOffset(id) << " " << layer << std::endl;
    // std::cout << "TC floats: "
    //           << "X/Z: " << tcFloat.getXOverZF() << ", "
    //           << "Y/Z: " << tcFloat.getYOverZF() << ", "
    //           << "Energy: " << tcFloat.getEnergyGeV() << std::endl;
      clusters_SA.push_back(tcFloat);
    // }
  }
}

void HGCalHistoClusteringWrapper::convertAlgorithmOutputs(std::vector<TPGCluster>& clusters_SA_out,
                                                          l1t::HGCalMulticlusterBxCollection& multiClusters_out,
                                                          const std::vector<edm::Ptr<l1t::HGCalCluster>>& tcPtrs) const {
  for (const auto& cluster : clusters_SA_out) {
    const auto hwCluster = cluster.getClData();
    auto sector = theConfiguration_.sector();
    double pt = cluster.getEnergyGeV();
    if ( pt < theConfiguration_.minClusterPtOut())
      continue;

    double phi = cluster.getGlobalPhiRad(sector);
    double eta = cluster.getGlobalEtaRad(sector);

    math::PtEtaPhiMLorentzVector clusterP4(pt, eta, phi, 0.);
    l1t::HGCalMulticluster multicluster;
    multicluster.setP4(clusterP4);
    multicluster.setMaxFinderPass(cluster.getMaxFinderPass());

    double sumLayerPt = 0;
    for ( const auto& itc : cluster.getCMSSWIndices() ) {
      const auto& tc_cmssw = tcPtrs.at(itc);
      sumLayerPt += tc_cmssw->pt();
      multicluster.addConstituent(tc_cmssw, false, 0.);
    }  
    //compute shower shapes
    shape_.fillShapes(multicluster, *triggerTools_.getTriggerGeometry());

    // Hardware cluster properties
    // Hardware value setters
    multicluster.setHwE(hwCluster.e);
    multicluster.setHwE_EM(hwCluster.e_em);
    multicluster.setHwGctBits(hwCluster.gctBits);
    multicluster.setHwFractionInCE_E(hwCluster.fractionInCE_E);
    multicluster.setHwFractionInCoreCE_E(hwCluster.fractionInCoreCE_E);
    multicluster.setHwFractionInEarlyCE_E(hwCluster.fractionInEarlyCE_E);
    multicluster.setHwFirstLayer(hwCluster.firstLayer);
    multicluster.setHwEta(hwCluster.w_eta);
    multicluster.setHwPhi(hwCluster.w_phi);
    multicluster.setHwZ(hwCluster.w_z);
    multicluster.setHwNTC(hwCluster.nTC);
    multicluster.setHwQualFlags(hwCluster.qualFlags);
    multicluster.setHwSigmaE(hwCluster.sigma_E);
    multicluster.setHwLastLayer(hwCluster.lastLayer);
    multicluster.setHwShowerLength(hwCluster.showerLength);
    multicluster.setHwSigmaZ(hwCluster.sigma_z);
    multicluster.setHwSigmaPhi(hwCluster.sigma_phi);
    multicluster.setHwCoreShowerLength(hwCluster.coreShowerLength);
    multicluster.setHwSigmaEta(hwCluster.sigma_eta);
    multicluster.setHwSigmaRoz(hwCluster.sigma_roz);

    // if ( pt > 10 ) {

    //   std::cout << "Got a cluster : " << pt << ", Eta : " << eta << ", Phi : " << phi
    //             << ", Sector : " << sector << ", Sum Layer pt: " << sumLayerPt
    //             << ", Constituent size: " << cluster.getCMSSWIndices().size() << std::endl;

    //   // Print first word members
    //   std::cout << "  e: " << hwCluster.e << ", e_em: " << hwCluster.e_em
    //         << ", gctBits: " << hwCluster.gctBits
    //         << ", fractionInCE_E: " << hwCluster.fractionInCE_E
    //         << ", fractionInCoreCE_E: " << hwCluster.fractionInCoreCE_E
    //         << ", fractionInEarlyCE_E: " << hwCluster.fractionInEarlyCE_E
    //         << ", firstLayer: " << hwCluster.firstLayer << std::endl;


    //   // Print second word members
    //   std::cout << "  w_eta: " << hwCluster.w_eta << ", w_phi: " << hwCluster.w_phi
    //         << ", w_z: " << hwCluster.w_z
    //         << ", nTC: " << hwCluster.nTC
    //         << ", qualFlags: " << hwCluster.qualFlags << std::endl;

    //   // Print third word members
    //   std::cout << "  sigma_E: " << hwCluster.sigma_E
    //         << ", lastLayer: " << hwCluster.lastLayer
    //         << ", showerLength: " << hwCluster.showerLength
    //         << ", sigma_z: " << hwCluster.sigma_z
    //         << ", sigma_phi: " << hwCluster.sigma_phi
    //         << ", coreShowerLength: " << hwCluster.coreShowerLength
    //         << ", sigma_eta: " << hwCluster.sigma_eta
    //         << ", sigma_roz: " << hwCluster.sigma_roz << std::endl;
    // }

    // double emIntfraction = l1thgcfirmware::Scales::floatFrac(hwCluster.fractionInCE_E);
    // multicluster.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::EM,
    //                                       emIntfraction * pt);

    // double emCoreIntfraction = l1thgcfirmware::Scales::floatFrac(hwCluster.fractionInCoreCE_E);
    // multicluster.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::EM_CORE,
    //                                       emCoreIntfraction * emIntfraction * pt);

    // double emHEarlyIntfraction = l1thgcfirmware::Scales::floatFrac(hwCluster.fractionInEarlyCE_E);
    // multicluster.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::H_EARLY,
    //                                       emHEarlyIntfraction * pt);


    // multicluster.setShowerLength(hwCluster.showerLength);
    // multicluster.setCoreShowerLength(hwCluster.coreShowerLength);
    // multicluster.setFirstLayer(hwCluster.firstLayer);
    // multicluster.setLast1layers(hwCluster.lastLayer);

    // multicluster.setZBarycenter(l1thgcfirmware::Scales::floatZ(hwCluster.w_z));
    // multicluster.setSigmaRRTot(l1thgcfirmware::Scales::floatSigmaRozRoz(hwCluster.sigma_roz));
    // multicluster.setSigmaEtaEtaTot(l1thgcfirmware::Scales::floatSigmaEta(hwCluster.sigma_eta));
    // multicluster.setSigmaPhiPhiTot(l1thgcfirmware::Scales::floatSigmaPhi(hwCluster.sigma_phi));
    // multicluster.setSigmaZZ(l1thgcfirmware::Scales::floatSigmaZ(hwCluster.sigma_z));
    // multicluster.setSigmaEE(l1thgcfirmware::Scales::floatSigmaE(hwCluster.sigma_E));

    // if ( pt > 90 && sumLayerPt < 80 ) {
    //   std::cout << "Cluster pt: " << pt << ", Sum Layer pt: " << sumLayerPt
    //             << ", Eta: " << eta << ", Phi: " << phi
    //             << ", Sector: " << sector
    //             << ", Constituent size: " << cluster.getCMSSWIndices().size() << std::endl;
    //   for ( const auto& itc : cluster.getCMSSWIndices() ) {
    //     const auto& tc_cmssw = tcPtrs.at(itc);
    //     std::cout << "  Constituent: " << tc_cmssw->pt() << ", Eta: " << tc_cmssw->eta() << ", Phi: " << tc_cmssw->phi()
    //               << ", Layer: " << triggerTools_.layerWithOffset(tc_cmssw->detId()) << std::endl;
    //   }
    // }
    // for (const auto& tc : cluster->constituents()) {
    //   const auto& tc_cmssw = inputClustersPtrs.at(tc->cmsswIndex().first).at(tc->cmsswIndex().second);
    //   // Add tc as constituent, but don't update any other properties of the multicluster i.e. leave them unchanged from those calculated by the emulator
    //   multicluster.addConstituent(tc_cmssw, false, 0.);
    // }


    multiClusters_out.push_back(0, multicluster);
    // if ( pt > 20 ) {
      // std::cout << "Got a cluster : " << pt << ", Eta : " << eta << ", Phi : " << phi << " " << sector << std::endl;
    // }
    // std::cout << "Got a cluster : " << pt << " " << eta << " " << phi << " " << sector << std::endl;
  }

  // for (const auto& cluster : clusterSums) {
  //   // Convert from digitised quantities
  //   if (cluster->w() == 0 || cluster->e() == 0)
  //     continue;
  //   double phi = (cluster->wphi() / cluster->w()) * theConfiguration_.phiRange() / theConfiguration_.phiNValues();
  //   double pt = cluster->e() / theConfiguration_.ptDigiFactor();

  //   if (pt < theConfiguration_.minClusterPtOut())
  //     continue;

  //   double rOverZ =
  //       (cluster->wroz() / cluster->w()) * theConfiguration_.rOverZRange() / theConfiguration_.rOverZNValues();
  //   double eta = -1.0 * std::log(tan(atan(rOverZ) / 2));
  //   eta *= theConfiguration_.zSide();

  //   auto sector = theConfiguration_.sector();
  //   phi = rotatePhiFromSectorZero(phi, sector);

  //   if (theConfiguration_.zSide() == 1) {
  //     phi = M_PI - phi;
  //   }
  //   phi -= (phi > M_PI) ? 2 * M_PI : 0;

  //   math::PtEtaPhiMLorentzVector clusterP4(pt, eta, phi, 0.);

  //   l1t::HGCalMulticluster multicluster;
  //   multicluster.setP4(clusterP4);

  //   // for (const auto& tc : cluster->constituents()) {
  //   //   const auto& tc_cmssw = inputClustersPtrs.at(tc->cmsswIndex().first).at(tc->cmsswIndex().second);
  //   //   // Add tc as constituent, but don't update any other properties of the multicluster i.e. leave them unchanged from those calculated by the emulator
  //   //   multicluster.addConstituent(tc_cmssw, false, 0.);
  //   // }

    // double emIntfraction = float(cluster->e_em()) / cluster->e();
    // multicluster.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::EM,
    //                                       emIntfraction * multicluster.energy());

    // double emCoreIntfraction = float(cluster->e_em_core()) / cluster->e();
    // multicluster.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::EM_CORE,
    //                                       emCoreIntfraction * multicluster.energy());

    // double emHEarlyIntfraction = float(cluster->e_h_early()) / cluster->e();
    // multicluster.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::H_EARLY,
    //                                       emHEarlyIntfraction * multicluster.energy());

  //   // Set cluster shower shape properties
  //   multicluster.setShowerLength(cluster->showerLen());
  //   multicluster.setCoreShowerLength(cluster->coreShowerLen());
  //   multicluster.setFirstLayer(cluster->firstLayer());
  //   multicluster.set_hw_sigma_e_quotient(cluster->sigma_e_quotient());
  //   multicluster.set_hw_sigma_e_fraction(cluster->sigma_e_fraction());
  //   multicluster.set_hw_mean_z_quotient(cluster->mean_z_quotient());
  //   multicluster.set_hw_mean_z_fraction(cluster->mean_z_fraction());
  //   multicluster.set_hw_mean_phi_quotient(cluster->mean_phi_quotient());
  //   multicluster.set_hw_mean_phi_fraction(cluster->mean_phi_fraction());
  //   multicluster.set_hw_mean_eta_quotient(cluster->mean_eta_quotient());
  //   multicluster.set_hw_mean_eta_fraction(cluster->mean_eta_fraction());
  //   multicluster.set_hw_mean_roz_quotient(cluster->mean_roz_quotient());
  //   multicluster.set_hw_mean_roz_fraction(cluster->mean_roz_fraction());
  //   multicluster.set_hw_sigma_z_quotient(cluster->sigma_z_quotient());
  //   multicluster.set_hw_sigma_z_fraction(cluster->sigma_z_fraction());
  //   multicluster.set_hw_sigma_phi_quotient(cluster->sigma_phi_quotient());
  //   multicluster.set_hw_sigma_phi_fraction(cluster->sigma_phi_fraction());
  //   multicluster.set_hw_sigma_eta_quotient(cluster->sigma_eta_quotient());
  //   multicluster.set_hw_sigma_eta_fraction(cluster->sigma_eta_fraction());
  //   multicluster.set_hw_sigma_roz_quotient(cluster->sigma_roz_quotient());
  //   multicluster.set_hw_sigma_roz_fraction(cluster->sigma_roz_fraction());
  //   multicluster.set_hw_e_em_over_e_quotient(cluster->e_em_over_e_quotient());
  //   multicluster.set_hw_e_em_over_e_fraction(cluster->e_em_over_e_fraction());
  //   multicluster.set_hw_e_em_core_over_e_em_quotient(cluster->e_em_core_over_e_em_quotient());
  //   multicluster.set_hw_e_em_core_over_e_em_fraction(cluster->e_em_core_over_e_em_fraction());
  //   multicluster.set_hw_e_h_early_over_e_quotient(cluster->e_h_early_over_e_quotient());
  //   multicluster.set_hw_e_h_early_over_e_fraction(cluster->e_h_early_over_e_fraction());

  //   multiClusters_out.push_back(0, multicluster);
  // }
}

void HGCalHistoClusteringWrapper::process(const std::vector<edm::Ptr<l1t::HGCalCluster>>& inputClusters,
                                          std::pair<l1t::HGCalMulticlusterBxCollection&, l1t::HGCalClusterBxCollection&>&
                                              outputMulticlustersAndRejectedClusters) const {
  std::vector<TPGTCBits> triggerCells_in_SA;
  convertCMSSWInputs(inputClusters, triggerCells_in_SA);
  std::vector<TPGCluster> clusters_SA_out;
  std::vector<l1thgcfirmware::HGCalCluster_HW> hwClusters_SA_out;
  makeClusters(triggerCells_in_SA, clusters_SA_out);
  convertAlgorithmOutputs(clusters_SA_out, outputMulticlustersAndRejectedClusters.first, inputClusters);
}

void HGCalHistoClusteringWrapper::makeClusters(const std::vector<TPGTCBits>& triggerCells_in_SA,
                                               std::vector<TPGCluster>& clusters_SA_out) const {

  TPGStage2Configuration::ClusPropLUT cplut;
  // cplut.readMuEtaLUT("/scratch/ec6821/HGC/Stage2Emu/ContractRush/CMSSW/CMSSW_14_2_2/src/L1Trigger/L1THGCal/data/mean_eta_LUT.csv");
  // cplut.readSigmaEtaLUT("/scratch/ec6821/HGC/Stage2Emu/ContractRush/CMSSW/CMSSW_14_2_2/src/L1Trigger/L1THGCal/data/sigma_eta_LUT.csv");
  cplut.readMuEtaLUT(muEtaLut_.fullPath().c_str());
  cplut.readSigmaEtaLUT(sigmaEtaLut_.fullPath().c_str());

  theAlgo_.setClusPropLUT(&cplut);
  // theAlgo_.setROverZ(theConfiguration_.sideLength() * sqrt(3.0));
  theAlgo_.run(triggerCells_in_SA, clusters_SA_out);
}

void HGCalHistoClusteringWrapper::configure(
    const std::tuple<const HGCalTriggerGeometryBase* const, const edm::ParameterSet&, const unsigned int, const int>&
        configuration) {

  setGeometry(std::get<0>(configuration));
  theConfiguration_.setSector(std::get<2>(configuration));
  theConfiguration_.setZSide(std::get<3>(configuration));

          // std::cout << "Configuring..." << std::endl;
    // TPGStage2Configuration::ClusPropLUT cplut;
    // cplut.readMuEtaLUT("/scratch/ec6821/HGC/Stage2Emu/ContractRush/CMSSW/CMSSW_14_2_2/src/L1Trigger/L1THGCal/data/mean_eta_LUT.csv");
    // cplut.readSigmaEtaLUT("/scratch/ec6821/HGC/Stage2Emu/ContractRush/CMSSW/CMSSW_14_2_2/src/L1Trigger/L1THGCal/data/sigma_eta_LUT.csv");
    // std::cout << "Read LUTs : " << cplut.muEtaSize() << " " << cplut.sigmaEtaSize() << std::endl;
    // theAlgo_.setClusPropLUT(&cplut);

          // setGeometry(std::get<0>(configuration));

  // theConfiguration_.setNTriggerLayers(std::get<0>(configuration)->lastTriggerLayer());
  // theConfiguration_.setTriggerLayers(std::get<0>(configuration)->triggerLayers());


  const edm::ParameterSet pset = std::get<1>(configuration)
                                     .getParameterSet("C3d_parameters")
                                     .getParameterSet("histoMax_C3d_clustering_parameters")
                                     .getParameterSet("layer2FwClusteringParameters");

  // Triangle size
  theConfiguration_.setSideLength(pset.getParameter<double>("sideLength"));
  // Parameters for selecting output clusters
  theConfiguration_.setMinClusterPtOut(pset.getParameter<double>("minClusterPtOut"));


  // theConfiguration_.setClusterizerOffset(pset.getParameter<unsigned int>("clusterizerOffset"));
  // theConfiguration_.setStepLatencies(pset.getParameter<std::vector<unsigned int>>("stepLatencies"));
  // theConfiguration_.setCClocks(pset.getParameter<unsigned int>("cClocks"));
  // theConfiguration_.setCInputs(pset.getParameter<unsigned int>("cInputs"));
  // theConfiguration_.setCInputs2(pset.getParameter<unsigned int>("cInputs2"));
  // theConfiguration_.setCInt(pset.getParameter<unsigned int>("cInt"));
  // theConfiguration_.setCColumns(pset.getParameter<unsigned int>("cColumns"));
  // theConfiguration_.setCRows(pset.getParameter<unsigned int>("cRows"));
  // theConfiguration_.setROverZHistOffset(pset.getParameter<unsigned int>("rOverZHistOffset"));
  // theConfiguration_.setROverZBinSize(pset.getParameter<unsigned int>("rOverZBinSize"));
  // theConfiguration_.setDepths(pset.getParameter<std::vector<unsigned int>>("depths"));
  // theConfiguration_.setLayerWeights_E(pset.getParameter<std::vector<unsigned int>>("layerWeights_E"));
  // theConfiguration_.setLayerWeights_E_EM(pset.getParameter<std::vector<unsigned int>>("layerWeights_E_EM"));
  // theConfiguration_.setLayerWeights_E_EM_core(pset.getParameter<std::vector<unsigned int>>("layerWeights_E_EM_core"));
  // theConfiguration_.setLayerWeights_E_H_early(pset.getParameter<std::vector<unsigned int>>("layerWeights_E_H_early"));
  // theConfiguration_.setCorrection(pset.getParameter<unsigned int>("correction"));
  // theConfiguration_.setSaturation(pset.getParameter<unsigned int>("saturation"));
  // const edm::ParameterSet& thresholdParams = pset.getParameterSet("thresholdMaximaParams");
  // theConfiguration_.setThresholdParams(thresholdParams.getParameter<unsigned int>("a"),
  //                                      thresholdParams.getParameter<unsigned int>("b"),
  //                                      thresholdParams.getParameter<int>("c"));

  // // Digitization parameters
  // const edm::ParameterSet& digitizationPset = pset.getParameterSet("digiParams");
  // theConfiguration_.setROverZRange(digitizationPset.getParameter<double>("rOverZRange"));
  // theConfiguration_.setROverZNValues(digitizationPset.getParameter<double>("rOverZNValues"));
  // theConfiguration_.setPhiRange(digitizationPset.getParameter<double>("phiRange"));
  // theConfiguration_.setPhiNValues(digitizationPset.getParameter<double>("phiNValues"));
  // theConfiguration_.setPtDigiFactor(digitizationPset.getParameter<double>("ptDigiFactor"));


  // // Input links parameters
  // const edm::ParameterSet& inputLinksPset = pset.getParameterSet("inputLinkParams");
  // theConfiguration_.setMaxClustersPerLink(inputLinksPset.getParameter<unsigned int>("maxClustersPerLink"));
  // theConfiguration_.setNInputLinks(inputLinksPset.getParameter<unsigned int>("nInputLinks"));

  // // TC distribution parameters
  // const edm::ParameterSet& tcDistPset = pset.getParameterSet("tcDistParams");
  // theConfiguration_.setN60Sectors(tcDistPset.getParameter<unsigned int>("n60Sectors"));
  // theConfiguration_.setNCoarsePhiDist1(tcDistPset.getParameter<unsigned int>("nCoarsePhiRegionsDist1"));
  // theConfiguration_.setNDistServers1(tcDistPset.getParameter<unsigned int>("nDistServers1"));
  // theConfiguration_.setDistServer1_nIn(tcDistPset.getParameter<unsigned int>("distServer1_nIn"));
  // theConfiguration_.setDistServer1_nOut(tcDistPset.getParameter<unsigned int>("distServer1_nOut"));
  // theConfiguration_.setDistServer1_nInterleave(tcDistPset.getParameter<unsigned int>("distServer1_nInterleave"));
  // theConfiguration_.setNCoarsePhiDist2(tcDistPset.getParameter<unsigned int>("nCoarsePhiRegionsDist2"));
  // theConfiguration_.setNDistServers2(tcDistPset.getParameter<unsigned int>("nDistServers2"));
  // theConfiguration_.setDistServer2_nIn(tcDistPset.getParameter<unsigned int>("distServer2_nIn"));
  // theConfiguration_.setDistServer2_nOut(tcDistPset.getParameter<unsigned int>("distServer2_nOut"));
  // theConfiguration_.setDistServer2_nInterleave(tcDistPset.getParameter<unsigned int>("distServer2_nInterleave"));

  // // Smearing parameters
  // const edm::ParameterSet& smearingPset = pset.getParameterSet("smearingParams");
  // theConfiguration_.setMaxBinsSmearing1D(smearingPset.getParameter<unsigned int>("maxBinsSmearing1D"));
  // theConfiguration_.setNBitsAreaNormLUT(smearingPset.getParameter<unsigned int>("nBitsAreaNormLUT"));

  // // Clusterizer parameters
  // const edm::ParameterSet& clusterizerPset = pset.getParameterSet("clusterizerParams");
  // theConfiguration_.setNBinsCosLUT(clusterizerPset.getParameter<unsigned int>("nBinsCosLUT"));
  // theConfiguration_.setNBitsCosLUT(clusterizerPset.getParameter<unsigned int>("nBitsCosLUT"));
  // theConfiguration_.setNFifos(clusterizerPset.getParameter<unsigned int>("nFifos"));
  // theConfiguration_.setNColumnsPerFifo(clusterizerPset.getParameter<unsigned int>("nColumnsPerFifo"));
  // theConfiguration_.setClusterizerMagicTime(clusterizerPset.getParameter<unsigned int>("clusterizerMagicTime"));
  // theConfiguration_.setFirstSeedBin(clusterizerPset.getParameter<unsigned int>("firstSeedBin"));
  // theConfiguration_.setNColumnFifoVeto(clusterizerPset.getParameter<unsigned int>("nColumnsFifoVeto"));
  // theConfiguration_.setDeltaR2Cut(clusterizerPset.getParameter<unsigned int>("deltaR2Cut"));
  // theConfiguration_.setNColumnsForClustering(clusterizerPset.getParameter<unsigned int>("nColumnsForClustering"));
  // theConfiguration_.setNRowsForClustering(clusterizerPset.getParameter<unsigned int>("nRowsForClustering"));

  // theConfiguration_.initializeLUTs();
};

double HGCalHistoClusteringWrapper::rotatePhiToSectorZero(const double phi, const unsigned sector) const {
  double rotatedPhi = phi;
  if (sector == 1) {
    if (rotatedPhi < M_PI and rotatedPhi > 0)
      rotatedPhi = rotatedPhi - (2. * M_PI / 3.);
    else
      rotatedPhi = rotatedPhi + (4. * M_PI / 3.);
  } else if (sector == 2) {
    rotatedPhi = rotatedPhi + (2. * M_PI / 3.);
  }
  return rotatedPhi;
}

double HGCalHistoClusteringWrapper::rotatePhiFromSectorZero(const double phi, const unsigned sector) const {
  double rotatedPhi = phi;
  if (sector == 1) {
    rotatedPhi += (2. * M_PI / 3.);
  } else if (sector == 2) {
    rotatedPhi += (4. * M_PI / 3.);
  }
  return rotatedPhi;
}

DEFINE_EDM_PLUGIN(HGCalHistoClusteringWrapperBaseFactory, HGCalHistoClusteringWrapper, "HGCalHistoClusteringWrapper");
