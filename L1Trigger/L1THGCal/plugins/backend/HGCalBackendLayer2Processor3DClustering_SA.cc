#include "L1Trigger/L1THGCal/interface/HGCalProcessorBase.h"

#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalHistoSeedingImpl.h"
#include "L1Trigger/L1THGCal/interface/HGCalAlgoWrapperBase.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalTriggerClusterInterpreterBase.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalStage2ClusterDistribution.h"
#include "L1Trigger/L1THGCal/interface/backend_semiemulator/TPGTCFloats.hh"

#include <utility>

class HGCalBackendLayer2Processor3DClusteringSA : public HGCalBackendLayer2ProcessorBase {
public:
  HGCalBackendLayer2Processor3DClusteringSA(const edm::ParameterSet& conf)
      : HGCalBackendLayer2ProcessorBase(conf),
        // distributor_(conf.getParameterSet("DistributionParameters")),
        conf_(conf) {


    // multiclusteringHistoSeeding_ = std::make_unique<HGCalHistoSeedingImpl>(
    //     conf.getParameterSet("C3d_parameters").getParameterSet("histoMax_C3d_seeding_parameters"));

    const edm::ParameterSet& clusteringParamConfig =
        conf.getParameterSet("C3d_parameters").getParameterSet("histoMax_C3d_clustering_parameters");
    const std::string& clusteringAlgoWrapperName = clusteringParamConfig.getParameter<std::string>("AlgoName");

    multiclusteringHistoClusteringWrapper_ = std::unique_ptr<HGCalHistoClusteringWrapperBase>{
        HGCalHistoClusteringWrapperBaseFactory::get()->create(clusteringAlgoWrapperName, clusteringParamConfig)};

    for (const auto& interpretationPset : conf.getParameter<std::vector<edm::ParameterSet>>("energy_interpretations")) {
      std::unique_ptr<HGCalTriggerClusterInterpreterBase> interpreter{
          HGCalTriggerClusterInterpreterFactory::get()->create(interpretationPset.getParameter<std::string>("type"))};
      interpreter->initialize(interpretationPset);
      energy_interpreters_.push_back(std::move(interpreter));
    }
  }

  void run(const edm::Handle<l1t::HGCalClusterBxCollection>& collHandle,
           std::pair<l1t::HGCalMulticlusterBxCollection, l1t::HGCalClusterBxCollection>& be_output) override {
    // if (multiclusteringHistoSeeding_) {
    //   multiclusteringHistoSeeding_->setGeometry(geometry());
    // }
    // if ( multiclusteringHistoClusteringWrapper_ ) {
    //   multiclusteringHistoClusteringWrapper_->setGeometry(geometry());
    // }

    l1t::HGCalMulticlusterBxCollection& collCluster3D_sorted = be_output.first;
    l1t::HGCalClusterBxCollection& rejectedClusters = be_output.second;

    /* create a persistent vector of pointers to the trigger-cells */
    std::unordered_map<uint32_t, std::vector<edm::Ptr<l1t::HGCalCluster>>> tcs_per_fpga;

    for (unsigned i = 0; i < collHandle->size(); ++i) {
      edm::Ptr<l1t::HGCalCluster> tc_ptr(collHandle, i);
      // if ( tc_ptr->position().z()< 0.0) {
      //   // Skip trigger cells with negative z position
      //   continue;
      // }
      // double eta = tc_ptr->eta();
      // double phi = tc_ptr->phi();
      // if (std::abs(eta) <= 2.5 || std::abs(eta) >= 2.7 || phi <= 0.3 || phi >= 0.6) {
      //   continue;
      // }
      for (uint32_t isect = 0 ; isect < 3 ; isect++ ){
        uint32_t addisect = 0;
        if(tc_ptr->position().z()>0.0) addisect = 3;
        double absZ = std::abs(tc_ptr->position().z());
        TPGTCFloats tcf0;
        tcf0.setROverZPhiF(tc_ptr->position().x()/absZ,tc_ptr->position().y()/absZ,isect+addisect);
        if(tcf0.getXOverZF()>=0.0){
          tcs_per_fpga[isect+addisect].push_back(tc_ptr);
          // std::cout << "[DEBUG] Added TC ptr with z, phi " << tc_ptr->position().z() << ", " << tc_ptr->phi() << " to tcs_per_fpga[" << isect+addisect << "]" << std::endl;

        }
      }

      // std::cout << "[DEBUG] Trigger cell ptr obtained" << std::endl;
      // uint32_t module = geometry()->getModuleFromTriggerCell(tc_ptr->detId());
      // std::cout << "[DEBUG] Module: " << module << std::endl;
      // uint32_t stage1_fpga = geometry()->getStage1FpgaFromModule(module);
      // std::cout << "[DEBUG] Stage1 FPGA: " << stage1_fpga << std::endl;
      // HGCalTriggerGeometryBase::geom_set stage2_fpgas = geometry()->getStage2FpgasFromStage1Fpga(stage1_fpga);

      // std::cout << "[DEBUG] stage2_fpgas size: " << stage2_fpgas.size() << std::endl;

      // uint16_t sec0(7), sec1(7);
      // auto phi_deg = tc_ptr->phi() * 180. / M_PI;
      // //The following segmentation is applied for trigger cells
      // if (phi_deg >= 0. and phi_deg <= 60.) {
      //   sec0 = 2;
      //   sec1 = 0;
      // } else if (phi_deg > 60. and phi_deg <= 120.) {
      //   sec0 = sec1 = 0;
      // } else if (phi_deg > 120. and phi_deg <= 180.) {
      //   sec0 = 0;
      //   sec1 = 1;
      // } else if (phi_deg >= -180. and phi_deg <= -120.) {
      //   sec0 = sec1 = 1;
      // } else if (phi_deg > -120. and phi_deg <= -60.) {
      //   sec0 = 1;
      //   sec1 = 2;
      // } else {
      //   sec0 = sec1 = 2;
      // }

      // if (tc_ptr->position().z() < 0.) {
      //   switch (sec0) {
      //     case 0:
      //       sec0 += 3;
      //       break;
      //     case 1:
      //       sec0 += 4;
      //       break;
      //     case 2:
      //       sec0 += 2;
      //       break;
      //     default:;
      //   }
      //   switch (sec1) {
      //     case 0:
      //       sec1 += 3;
      //       break;
      //     case 1:
      //       sec1 += 4;
      //       break;
      //     case 2:
      //       sec1 += 2;
      //       break;
      //     default:;
      //   }
      // }

      // unsigned int addSec = (tc_ptr->position().z() > 0) ? 0 : 3;
      // tcs_per_fpga[sec0 + addSec].push_back(tc_ptr);
      // if (sec1 != sec0) {
      //   tcs_per_fpga[sec1 + addSec].push_back(tc_ptr);
      // }
      // std::cout << "[DEBUG] Added cluster ptr with phi " << tc_ptr->phi() << " to tcs_per_fpga[" << sec0 << " " << sec0 + addSec << " "<< sec1 << " " << sec1 + addSec << "]" << std::endl;
      // for (auto& fpga : stage2_fpgas) {
      //   tcs_per_fpga[fpga].push_back(tc_ptr);
      //   std::cout << "[DEBUG] Added cluster ptr to tcs_per_fpga[" << fpga << "]" << std::endl;
      // }
    }

    // Configuration
    const std::pair<const HGCalTriggerGeometryBase* const, const edm::ParameterSet&> configuration{geometry(), conf_};

    for (auto& fpga_tcs : tcs_per_fpga) {
      // Inputs
      const std::vector<edm::Ptr<l1t::HGCalCluster>>& inputClusters_perFPGA{fpga_tcs.second};
      // Outputs
      l1t::HGCalMulticlusterBxCollection collCluster3D_perFPGA;
      l1t::HGCalClusterBxCollection rejectedClusters_perFPGA;

      std::pair<l1t::HGCalMulticlusterBxCollection&, l1t::HGCalClusterBxCollection&>
          outputMulticlustersAndRejectedClusters_perFPGA{collCluster3D_perFPGA, rejectedClusters_perFPGA};

      HGCalTriggerBackendDetId stage2_fpga_id(fpga_tcs.first);
      // const auto stage2_sector = stage2_fpga_id.sector();
      const auto stage2_sector = fpga_tcs.first;
      const auto zSide = stage2_fpga_id.zside();
      const auto clusteringConfig = std::make_tuple(geometry(), std::ref(conf_), stage2_sector, zSide);



      if (multiclusteringHistoClusteringWrapper_) {
        // std::cout << "Configuring clustering wrapper for FPGA " << stage2_fpga_id << std::endl;
        multiclusteringHistoClusteringWrapper_->configure(clusteringConfig);
      }

      // Process
      multiclusteringHistoClusteringWrapper_->process(inputClusters_perFPGA,
                                                      outputMulticlustersAndRejectedClusters_perFPGA);

      for (const auto& collcluster : collCluster3D_perFPGA) {
        collCluster3D_sorted.push_back(0, collcluster);
      }
      for (const auto& rejectedcluster : rejectedClusters_perFPGA) {
        rejectedClusters.push_back(0, rejectedcluster);
      }
    }
  }

private:
  /* algorithms instances */
  // std::unique_ptr<HGCalHistoSeedingImpl> multiclusteringHistoSeeding_;

  std::unique_ptr<HGCalHistoClusteringWrapperBase> multiclusteringHistoClusteringWrapper_;

  std::vector<std::unique_ptr<HGCalTriggerClusterInterpreterBase>> energy_interpreters_;

  // HGCalStage2ClusterDistribution distributor_;
  const edm::ParameterSet conf_;

};

DEFINE_EDM_PLUGIN(HGCalBackendLayer2Factory,
                  HGCalBackendLayer2Processor3DClusteringSA,
                  "HGCalBackendLayer2Processor3DClusteringSA");
