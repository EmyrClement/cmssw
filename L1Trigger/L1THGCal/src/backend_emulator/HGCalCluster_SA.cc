#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"
#include <algorithm>
#include <cmath>

// Many useful functions for handling and manipulation of bitsets
// From L1 Track Finder DTC
// But could perhaps still just use ap_uint?
#include "DataFormats/L1TrackTrigger/interface/TTBV.h"

using namespace l1thgcfirmware;

const HGCalCluster& HGCalCluster::operator+=(HGCalCluster& c) {
  // Not handling field widths
  HGCalCluster original(*this);
  this->set_n_tc(this->n_tc() + c.n_tc());
  this->set_e(this->e() + c.e());
  this->set_e_em(this->e_em() + c.e_em());
  this->set_e_em_core(this->e_em_core() + c.e_em_core());
  this->set_e_h_early(this->e_h_early() + c.e_h_early());
  this->set_w(this->w() + c.w());
  this->set_n_tc_w(this->n_tc_w() + c.n_tc_w());
  this->set_w2(this->w2() + c.w2());
  this->set_wz(this->wz() + c.wz());
  this->set_weta(this->weta() + c.weta());
  this->set_wphi(this->wphi() + c.wphi());
  this->set_wroz(this->wroz() + c.wroz());
  this->set_wz2(this->wz2() + c.wz2());
  this->set_weta2(this->weta2() + c.weta2());
  this->set_wphi2(this->wphi2() + c.wphi2());
  this->set_wroz2(this->wroz2() + c.wroz2());

  this->set_layerbits(this->layerbits() | c.layerbits());
  this->set_sat_tc(this->sat_tc() | c.sat_tc());
  this->set_shapeq(this->shapeq() | c.shapeq());

  const unsigned clusterWeightSat80 = ((1 << 16) - 1) * 0.8;  // 52428
  if (w_ <= clusterWeightSat80 && original.shapeq() == 1 && c.shapeq() == 1) {
    this->set_shapeq(1);
  } else {
    this->set_shapeq(0);

    if (this->w() > c.w()) {
      this->set_w(original.w());
      this->set_w2(original.w2());
      this->set_wz(original.wz());
      this->set_weta(original.weta());
      this->set_wphi(original.wphi());
      this->set_wroz(original.wroz());
      this->set_wz2(original.wz2());
      this->set_weta2(original.weta2());
      this->set_wphi2(original.wphi2());
      this->set_wroz2(original.wroz2());
      this->set_n_tc_w(original.n_tc_w());
    } else {
      this->set_w(c.w());
      this->set_w2(c.w2());
      this->set_wz(c.wz());
      this->set_weta(c.weta());
      this->set_wphi(c.wphi());
      this->set_wroz(c.wroz());
      this->set_wz2(c.wz2());
      this->set_weta2(c.weta2());
      this->set_wphi2(c.wphi2());
      this->set_wroz2(c.wroz2());
      this->set_n_tc_w(c.n_tc_w());
    }
  }

  for (const auto& constituent : c.constituents()) {
    this->add_constituent(constituent);
  }

  return *this;
}

void HGCalCluster::clearClusterWords() {
  for ( auto& clusterWord : packedData_ ) {
    clusterWord = 0;
  }
}

HGCalCluster::ClusterWords HGCalCluster::formatClusterWords( const ClusterAlgoConfig& config ) {
  clearClusterWords();
  packedData_[0] = formatFirstWord( config );
  packedData_[1] = formatSecondWord( config );
  packedData_[2] = formatThirdWord( config );
  packedData_[3] = formatFourthWord( config );

  return packedData_;
}

HGCalCluster::ClusterWord HGCalCluster::formatFirstWord( const ClusterAlgoConfig& config ) {

  const unsigned nb_e = 14;
  const unsigned nb_e_em = 14;
  const unsigned nb_gctEGSelectBits = 4;
  const unsigned nb_fractionInCE_E = 8;
  const unsigned nb_fractionInCoreCE_E = 8;
  const unsigned nb_fractionInEarlyCE_H = 8;
  const unsigned nb_firstLayer = 6;
  const unsigned nb_spare = 2;

  unsigned et = round(float(e())/256);
  ap_uint<nb_e> hw_e = et;
  unsigned et_em = round(float(e_em())/256);
  ap_uint<nb_e_em> hw_e_em = et_em;

  ap_uint<nb_fractionInCE_E> hw_fractionInCE_E = round( 128 * float(e_em()) / e() );
  ap_uint<nb_fractionInCoreCE_E> hw_fractionInCoreCE_E = round( 128 * float(e_em_core()) / e_em() );
  ap_uint<nb_fractionInEarlyCE_H> hw_fractionInEarlyCE_H = round( 128 * float(e_h_early()) / e() );

  ap_uint<1> gctBit0 = hw_fractionInCE_E > 64;
  ap_uint<1> gctBit1 = hw_fractionInCoreCE_E > 64;
  ap_uint<1> gctBit2 = hw_fractionInEarlyCE_H > 64;
  ap_uint<1> gctBit3 = et_em > 64;
  ap_uint<4> gctSelectBits = (gctBit3, gctBit2, gctBit1, gctBit0);

  ap_uint<nb_firstLayer> hw_firstLayer = firstLayer();
  ap_uint<nb_spare> hw_spare = 0;


  const ap_uint<wordLength> clusterWord = (
    hw_spare,
    hw_firstLayer,
    hw_fractionInEarlyCE_H,
    hw_fractionInCoreCE_E,
    hw_fractionInCE_E,
    gctSelectBits,
    hw_e_em,
    hw_e
  );

  return clusterWord.to_ulong();
}
HGCalCluster::ClusterWord HGCalCluster::formatSecondWord( const ClusterAlgoConfig& config ) {

  const unsigned nb_eta = 9;
  const unsigned nb_phi = 9;
  const unsigned nb_z = 14;
  const unsigned nb_nTC = 10;
  const unsigned nb_satTC = 1;
  const unsigned nb_qualFracCE_E = 1;
  const unsigned nb_qualFracCoreCE_E = 1;
  const unsigned nb_qualFracEarlyCE_H = 1;
  const unsigned nb_qualSigmasMeans = 1;
  const unsigned nb_satPhi = 1;
  const unsigned nb_nominalPhi = 1;
  const unsigned nb_spare = 15;

  ap_uint<nb_eta> hw_eta = 317;// round( float(weta()) / w() ) + 317;
  ap_int<nb_phi> hw_phi = 0;
  ap_uint<nb_z> hw_z = 0;
  ap_uint<nb_nTC> hw_nTC = n_tc();
  ap_uint<nb_satTC> hw_satTC = sat_tc();
  ap_uint<nb_qualFracCE_E> hw_qualFracCE_E = e() != 0x3FFFFF && e_em() != 0x3FFFFF;
  ap_uint<nb_qualFracCoreCE_E> hw_qualFracCoreCE_E = e_em_core() != 0x3FFFFF && e_em() != 0x3FFFFF;
  ap_uint<nb_qualFracEarlyCE_H> hw_qualFracEarlyCE_H = e_h_early() != 0x3FFFFF && e() != 0x3FFFFF;
  ap_uint<nb_qualSigmasMeans> hw_qualSigmasMeans = shapeq();
  ap_uint<nb_satPhi> hw_satPhi = 0;
  ap_uint<nb_nominalPhi> hw_nominalPhi = 0;
  ap_uint<nb_spare> hw_spare = 0;

  const int hw_phi10b = int( (float(wphi()) / w()) * (5./24) ) - 360;
  if( hw_phi10b > 255 ) {
    hw_phi = 255;
    hw_satPhi = 1;
  }
  else if ( hw_phi10b < -256 ) {
    hw_phi = -256;
    hw_satPhi = 1;
  }
  else {
    hw_phi = hw_phi10b;
    hw_satPhi = 0;
  }
  hw_nominalPhi = ( hw_phi10b < 240 ) && ( hw_phi10b > -241 );
  std::cout << "E, nTC : " << e() << " " << n_tc() << std::endl;
  std::cout << "Eta : " << hw_eta << " " << hw_eta.to_string() << std::endl;
  std::cout << "Phi : " << hw_phi << " " << hw_phi.to_string() << std::endl;
  std::cout << "Sat phi : " << hw_satPhi << std::endl;
  std::cout << "Nominal phi : " << hw_nominalPhi << std::endl;


  hw_z = round( float(wz()) / w() ) + ( 3221 * 2 );

  const ap_uint<wordLength> clusterWord = (
    hw_spare,
    hw_nominalPhi,
    hw_satPhi,
    hw_qualSigmasMeans,
    hw_qualFracEarlyCE_H,
    hw_qualFracCoreCE_E,
    hw_qualFracCE_E,
    hw_satTC,
    hw_nTC,
    hw_z,
    hw_phi,
    hw_eta
  );
  std::cout << "Cluster word : " << clusterWord << " " << clusterWord.to_string() << std::endl;
  // // Not being careful with round/floor here
  // constexpr float ETAPHI_LSB = M_PI / 720;
  // double rOverZ = ( 1.0 * wroz() / w() ) * config.rOverZRange() / config.rOverZNValues();
  // double eta = -1.0 * std::log( tan( atan( rOverZ ) / 2 ) );
  // // eta *= theConfiguration_.zSide();

  // // Interface document definition (old?)
  // // eta -= 1.5; 
  // // int l1Eta = eta / ETAPHI_LSB;
  // // const TTBV hw_Eta(l1Eta, 9, false);
  // // double etaTemp = eta - 1.5;
  // // const TTBV hw_Eta_interface(int(etaTemp / ETAPHI_LSB), 9, false);

  // // Current correlator definition
  // // eta -= 2.25;
  // // int l1Eta = round( eta / ETAPHI_LSB );
  // // const TTBV hw_Eta(l1Eta, 9, true);
  // // std::cout << "Global eta, l1 eta : " << -1.0 * std::log( tan( atan( rOverZ ) / 2 ) ) << " " << l1Eta << " " << hw_Eta.str() << std::endl;
  // // int l1Phi = ( ( 1.0 * wphi() / w() ) * config.phiRange() / config.phiNValues() - M_PI/3 )  / ETAPHI_LSB;

  // const TTBV hw_Eta(0, 9, true);


  // // double phi = ( 1.0 * wphi() / w() ) * config.phiRange() / config.phiNValues();
  // // if ( config.zSide() == 1 ) {
  // //   phi = M_PI - phi;
  // // }
  // // // phi -= ( phi > M_PI ) ? 2 * M_PI : 0;

  // // int l1Phi = round( ( phi - M_PI/2  )  / ETAPHI_LSB );
  // // const TTBV hw_Phi(l1Phi, 9, true);

  // const TTBV hw_Phi(0, 9, true);
  // const TTBV hw_z(0, 12);
  // const TTBV hw_spare1(0, 2);

  // const TTBV hw_nTCs(int(n_tc()), 10);
  // // const TTBV hw_nTCs(0, 10);
  // const TTBV hw_qualityFlags(0, 10);
  // const TTBV hw_spare2(0, 12);

  // // ClusterWord clusterWord = TTBV(hw_Eta + hw_Phi + hw_z + hw_spare1 + hw_nTCs + hw_qualityFlags + hw_spare2).bs();
  // ClusterWord clusterWord = TTBV(    
  //   hw_spare2 +
  //   hw_qualityFlags +
  //   hw_nTCs + 
  //   hw_spare1 + 
  //   hw_z + 
  //   hw_Phi + 
  //   hw_Eta
  //   ).bs();

  return clusterWord.to_ulong();
}
HGCalCluster::ClusterWord HGCalCluster::formatThirdWord( const ClusterAlgoConfig& config ) {
  ClusterWord clusterWord = 0;
  return clusterWord;
}

HGCalCluster::ClusterWord HGCalCluster::formatFourthWord( const ClusterAlgoConfig& config ) {
  ClusterWord clusterWord = 0;
  return clusterWord;
}

void HGCalCluster::clearClusterSumWords() {
  for ( auto& clusterSumWord : packedData_clustersSums_ ) {
    clusterSumWord = 0;
  }
}

HGCalCluster::ClusterSumWords HGCalCluster::formatClusterSumWords( const ClusterAlgoConfig& config ) {
  clearClusterSumWords();

  const unsigned nb_nTCs = 10;
  const unsigned nb_e = 22;
  const unsigned nb_e_em = 22;
  const unsigned nb_e_em_core = 22;
  const unsigned nb_e_h_early = 22;

  const unsigned nb_w = 16;
  const unsigned nb_n_tc_w = 10;
  const unsigned nb_w2 = 32;
  const unsigned nb_wz = 29;
  const unsigned nb_weta = 26;
  const unsigned nb_wphi = 28;
  const unsigned nb_wroz = 29;

  const unsigned nb_wz2 = 42;
  const unsigned nb_weta2 = 36;
  const unsigned nb_wphi2 = 40;
  const unsigned nb_wroz2 = 42;
  const unsigned nb_layerbits = 36;
  const unsigned nb_sat_tc = 1;
  const unsigned nb_shapeq = 1;

  ap_uint<nb_nTCs> hw_nTCs = n_tc();
  ap_uint<nb_e> hw_e = e();
  ap_uint<nb_e_em> hw_e_em = e_em();
  ap_uint<nb_e_em_core> hw_e_em_core = e_em_core();
  ap_uint<nb_e_h_early> hw_e_h_early = e_h_early();

  ap_uint<nb_w> hw_w = w();
  ap_uint<nb_n_tc_w> hw_n_tc_w = n_tc_w();
  ap_uint<nb_w2> hw_w2 = 0; //w2();
  ap_uint<nb_wz> hw_wz = 0;
  ap_uint<nb_weta> hw_weta = 0;
  ap_uint<nb_wphi> hw_wphi = wphi();
  ap_uint<nb_wroz> hw_wroz = 0;

  ap_uint<nb_wz2> hw_wz2 = 0;
  ap_uint<nb_weta2> hw_weta2 = 0;
  ap_uint<nb_wphi2> hw_wphi2 = 0;
  ap_uint<nb_wroz2> hw_wroz2 = 0;
  ap_uint<nb_layerbits> hw_layerbits = 0;
  ap_uint<nb_sat_tc> hw_sat_tc = 0;
  ap_uint<nb_shapeq> hw_shapeq = 0;

  const ap_uint<allClusterSumWordsLength> clusterSumRecord = (
    hw_shapeq,
    hw_sat_tc,
    hw_layerbits,
    hw_wroz2,
    hw_wphi2,
    hw_weta2,
    hw_wz2,
    hw_wroz,
    hw_wphi,
    hw_weta,
    hw_wz,
    hw_w2,
    hw_n_tc_w,
    hw_w,
    hw_e_h_early,
    hw_e_em_core,
    hw_e_em,
    hw_e,
    hw_nTCs
  );

  for ( unsigned iWord = 0; iWord < nWordsPerClusterSum; ++iWord ) {
    ClusterSumWord word = clusterSumRecord.range((iWord+1)*(clusterSumWordLength)-1,iWord*clusterSumWordLength).to_ulong();
    packedData_clustersSums_[iWord] = word;
  }
  // std::cout << "=== Cluster info ===" << std::endl;
  // std::cout << "NTCs : " << n_tc() << " " << hw_nTCs << " " << hw_nTCs.to_string() << std::endl;
  // std::cout << "E : " << e() << " " << hw_e << " " << hw_e.to_string() << std::endl;
  // std::cout << "w : " << w() << " " << hw_w << " " << hw_w.to_string() << std::endl;
  // std::cout << "n_tc_w : " << n_tc_w() << " " << hw_n_tc_w << " " << hw_n_tc_w.to_string() << std::endl;
  // std::cout << "w2 : " << w2() << " " << hw_w2 << " " << hw_w2.to_string() << std::endl;
  // std::cout << "Cluster sum record : " << clusterSumRecord << " " << clusterSumRecord.to_string()<< std::endl;

  return packedData_clustersSums_;
}