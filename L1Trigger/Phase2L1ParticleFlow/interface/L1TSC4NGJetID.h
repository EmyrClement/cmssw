#ifndef L1TRIGGER_PHASE2L1PARTICLEFLOWS_L1TSC4NGJetID_H
#define L1TRIGGER_PHASE2L1PARTICLEFLOWS_L1TSC4NGJetID_H

#include <string>
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/PFJet.h"
#include "DataFormats/L1TParticleFlow/interface/datatypes.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/jetmet/L1SeedConePFJetEmulator.h"

//HLS4ML compiled emulator modeling
#include "ap_fixed.h"
#include "hls4ml/emulator.h"

namespace L1TSC4NGJet{

  constexpr int ceillog2(int x){
    return (x <= 2) ? 1 : 1 + ceillog2((x+1) / 2);
  }

  template<class data_T, int N>
  inline float real_val_from_idx(unsigned i){
      // Treat the index as the top N bits
      static constexpr int NB = ceillog2(N); // number of address bits for table
      data_T x(0);
      // The MSB of 1 is implicit in the table
      x[x.width-1] = 1;
      // So we can use the next NB bits for real data
      x(x.width-2, x.width-NB-1) = i;
      return (float) x;
  }

  template<class data_T, int N>
  inline unsigned idx_from_real_val(data_T x){
      // Slice the top N bits to get an index into the table
      static constexpr int NB = ceillog2(N); // number of address bits for table
      // Slice the top-1 NB bits of the value
      // the MSB of '1' is implicit, so only slice below that
      ap_uint<NB> y = x(x.width-2, x.width-NB-1);
      return (unsigned) y(NB-1, 0);
  }


  template<class data_T, class table_T, int N>
  void init_invert_table(table_T table_out[N]){
    // The template data_T is the data type used to address the table
    for(unsigned i = 0; i < N; i++){
        float x = real_val_from_idx<data_T, N>(i);
        table_T inv_x = 1 / x;
        table_out[i] = inv_x;
    }
  }

  template<class in_t, class table_t, int N>
  table_t invert_with_shift(in_t in){
    table_t inv_table[N];
    init_invert_table<in_t, table_t, N>(inv_table);

    // find the first '1' in the denominator
    int msb = 0;
    for(int b = 0; b < in.width; b++){
        if(in[b]) msb = b;
    }
    // shift up the denominator such that the left-most bit (msb) is '1'
    in_t in_shifted = in << (in.width-msb-1);
    // lookup the inverse of the shifted input
    int idx = idx_from_real_val<in_t,N>(in_shifted);
    table_t inv_in = inv_table[idx];
    // shift the output back
    table_t out = inv_in << (in.width-msb-1);
    return out;
}

template<class t>
t candidate_mass(l1ct::PuppiObj puppicand) {
  // Define lookup table
  static const t PION_MASS = t(0.13);
  static const t PHOTON_MASS = t(0.0);
  static const t ELECTRON_MASS = t(0.005);
  static const t MUON_MASS = t(0.105);
  static const t K_MASS = t(0.5);

  // Default to pion mass
  t massCand = PION_MASS;

  if (puppicand.hwId.bits == l1ct::ParticleID::PHOTON) {
    massCand = PHOTON_MASS;
  }
  else if (puppicand.hwId.bits == l1ct::ParticleID::ELEPLUS || 
           puppicand.hwId.bits == l1ct::ParticleID::ELEMINUS) {
    massCand = ELECTRON_MASS;
  }
  else if (puppicand.hwId.bits == l1ct::ParticleID::MUMINUS || 
           puppicand.hwId.bits == l1ct::ParticleID::MUPLUS) {
    massCand = MUON_MASS;
  }
  else if (puppicand.hwId.bits == l1ct::ParticleID::HADZERO) {
    massCand = K_MASS;
  }

  return massCand;
}

}

class L1TSC4NGJetID {
public:
  L1TSC4NGJetID(const std::shared_ptr<hls4mlEmulator::Model> model, int iNParticles, bool debug);

  typedef ap_fixed<24, 12, AP_RND, AP_SAT, 0> inputtype;
  typedef std::array<ap_ufixed<24, 12, AP_RND, AP_SAT, 0>, 8> classtype;
  typedef std::array<ap_fixed<16, 6>, 1> regressiontype;
  typedef std::pair<regressiontype, classtype> pairtype;

  void setNNVectorVar();
  std::vector<float> EvaluateNNFixed();
  std::vector<float> computeFixed(const l1t::PFJet &iJet, bool useRawPt);

private:
  std::vector<inputtype> NNvectorVar_;
  int fNParticles_;
  unique_ptr<inputtype[]> fPt_;
  unique_ptr<inputtype[]> fPt_rel_;
  unique_ptr<inputtype[]> fDEta_;
  unique_ptr<inputtype[]> fDPhi_;
  unique_ptr<inputtype[]> fPt_log_;
  unique_ptr<inputtype[]> fMass_;
  unique_ptr<inputtype[]> fZ0_;
  unique_ptr<inputtype[]> fDxy_;
  unique_ptr<inputtype[]> fIs_filled_;
  unique_ptr<inputtype[]> fPuppi_weight_;
  unique_ptr<inputtype[]> fEmID_;
  unique_ptr<inputtype[]> fQuality_;

  unique_ptr<inputtype[]> fCharge_;
  unique_ptr<inputtype[]> fId_;
  std::shared_ptr<hls4mlEmulator::Model> modelRef_;

  bool isDebugEnabled_;
};
#endif
