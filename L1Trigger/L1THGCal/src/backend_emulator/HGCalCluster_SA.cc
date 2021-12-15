#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"
#include <math.h>
#include <algorithm>


using namespace l1thgcfirmware;

double HGCalCluster::Sigma_Energy(unsigned int N_TC_W, unsigned long int Sum_W2, unsigned int Sum_W) {
  unsigned long int N = N_TC_W*Sum_W2 - pow(Sum_W,2);
  unsigned long int D = pow(N_TC_W,2);
  return  (double)sqrt(N/D);
}

double HGCalCluster::Mean_coordinate(unsigned int Sum_Wc, unsigned int Sum_W) {
  return (double)Sum_Wc/Sum_W;
}

double HGCalCluster::Sigma_Coordinate(unsigned int Sum_W, unsigned long int Sum_Wc2, unsigned int Sum_Wc) {
  unsigned long int N = Sum_W*Sum_Wc2 - pow(Sum_Wc,2);
  unsigned long int D = pow(Sum_W,2);
  return (double)sqrt(N/D);
}

double HGCalCluster::Energy_ratio(unsigned int E_N, unsigned int E_D) {
  if ( E_D == 0 ) { 
    return 0; 
  } else {
    return (double)E_N/E_D;
  }
}

std::vector<int> HGCalCluster::ShowerLengthProperties(unsigned long int layerBits) {
  int counter = 0;
  int firstLayer = 0;
  bool firstLayerFound = false;
  int lastLayer = 0;
 

  std::vector<int> layerBits_array;
  for (int idx = 0; idx<36; idx++) {
    if ( (layerBits&(1L<<(35-idx)) ) >= 1L ) {
      if ( !firstLayerFound ) {
        firstLayer = idx+1;
        firstLayerFound=true;
      }
      lastLayer = idx+1;
      counter += 1;
    } else {
      layerBits_array.push_back(counter);
      counter = 0;
    }
  }
  int showerLen = lastLayer - firstLayer + 1;
  int coreShowerLen = 36;
  if ( layerBits_array.size()>0 ) {
    coreShowerLen = *std::max_element(layerBits_array.begin(), layerBits_array.end()); 
  }

  std::vector<int> output = {firstLayer,lastLayer, showerLen, coreShowerLen};

  return output;
}


unsigned long int HGCalCluster::Sigma_E_Fraction()  { 
  double Sigma_E_temp = Sigma_Energy( this->n_tc_w(), this->w2(), this->w() );
  double intpart;
  return modf(Sigma_E_temp,&intpart)*pow(2,1);
}


unsigned long int HGCalCluster::Sigma_E_Quotient()  { 
  double Sigma_E_temp = Sigma_Energy( this->n_tc_w(), this->w2(), this->w() );
  double intpart;
  double frac =  modf(Sigma_E_temp,&intpart);
  return intpart;
}

unsigned long int HGCalCluster::Mean_z_Fraction()   { 
  double Mean_z_temp = Mean_coordinate(this->wz(), this->w());
  double intpart;
  return modf(Mean_z_temp,&intpart)*pow(2,2);
}

unsigned long int HGCalCluster::Mean_z_Quotient()  { 
  double Mean_z_temp = Mean_coordinate(this->wz(), this->w());
  double intpart;
  double frac =  modf(Mean_z_temp,&intpart)*pow(2,2);
  return intpart;
}
    
unsigned long int HGCalCluster::Mean_phi_Fraction()  { 
  double Mean_phi_temp = Mean_coordinate(this->wphi(), this->w());
  double intpart;
  return  modf(Mean_phi_temp,&intpart)*pow(2,2);
}

unsigned long int HGCalCluster::Mean_phi_Quotient()  { 
  double Mean_phi_temp = Mean_coordinate(this->wphi(), this->w());
  double intpart;
  double frac =  modf(Mean_phi_temp,&intpart)*pow(2,2);
  return intpart;
}
    
unsigned long int HGCalCluster::Mean_eta_Fraction() { 
  double Mean_eta_temp = Mean_coordinate(this->weta(), this->w());
  double intpart;
  return  modf(Mean_eta_temp,&intpart)*pow(2,2);
}
   
unsigned long int HGCalCluster::Mean_eta_Quotient()  { 
  double Mean_eta_temp = Mean_coordinate(this->weta(), this->w());
  double intpart;
  double frac =  modf(Mean_eta_temp,&intpart)*pow(2,2);
  return intpart;
}
    
unsigned long int HGCalCluster::Mean_roz_Fraction() { 
  double Mean_roz_temp = Mean_coordinate(this->wroz(), this->w());
  double intpart;
  return  modf(Mean_roz_temp,&intpart)*pow(2,2);
}

unsigned long int HGCalCluster::Mean_roz_Quotient()  { 
  double Mean_roz_temp = Mean_coordinate(this->wroz(), this->w());
  double intpart;
  double frac =  modf(Mean_roz_temp,&intpart)*pow(2,2);
  return intpart;
}

unsigned long int HGCalCluster::Sigma_z_Fraction()  { 
  double Sigma_z_temp = Sigma_Coordinate(this->w(), this->wz2(), this->wz());
  double intpart;  
  return modf(Sigma_z_temp,&intpart)*pow(2,1);
}

unsigned long int HGCalCluster::Sigma_z_Quotient()  { 
  double Sigma_z_temp = Sigma_Coordinate(this->w(), this->wz2(), this->wz());
  double intpart;
  double frac = modf(Sigma_z_temp,&intpart)*pow(2,1);
  return intpart;
}

unsigned long int HGCalCluster::Sigma_phi_Fraction()  { 
  double Sigma_phi_temp = Sigma_Coordinate(this->w(), this->wphi2(), this->wphi());
  double intpart;
  return modf(Sigma_phi_temp,&intpart)*pow(2,1);
}

unsigned long int HGCalCluster::Sigma_phi_Quotient()  { 
  double Sigma_phi_temp = Sigma_Coordinate(this->w(), this->wphi2(), this->wphi());
  double intpart;
  double frac = modf(Sigma_phi_temp,&intpart)*pow(2,1);
  return intpart;
}

unsigned long int HGCalCluster::Sigma_eta_Fraction()  { 
  double Sigma_eta_temp = Sigma_Coordinate(this->w(), this->weta2(), this->weta());
  double intpart;
  return modf(Sigma_eta_temp,&intpart)*pow(2,1);
}

unsigned long int HGCalCluster::Sigma_eta_Quotient()  { 
  double Sigma_eta_temp = Sigma_Coordinate(this->w(), this->weta2(), this->weta());
  double intpart;
  double frac = modf(Sigma_eta_temp,&intpart)*pow(2,1);
  return intpart;
}

unsigned long int HGCalCluster::Sigma_roz_Fraction()  { 
  double Sigma_roz_temp = Sigma_Coordinate(this->w(), this->wroz2(), this->wroz());
  double intpart;
  return modf(Sigma_roz_temp,&intpart)*pow(2,1);
}

unsigned long int HGCalCluster::Sigma_roz_Quotient()  { 
  double Sigma_roz_temp = Sigma_Coordinate(this->w(), this->wroz2(), this->wroz());
  double intpart;
  double frac = modf(Sigma_roz_temp,&intpart)*pow(2,1);
  return intpart;
}

unsigned long int HGCalCluster::FirstLayer()  { 
  std::vector<int> layeroutput = ShowerLengthProperties(this->layerbits());
  return layeroutput[0];
}

unsigned long int HGCalCluster::LastLayer()  { 
  std::vector<int> layeroutput = ShowerLengthProperties(this->layerbits());
  return layeroutput[1];
}

unsigned long int HGCalCluster::ShowerLen()  { 
  std::vector<int> layeroutput = ShowerLengthProperties(this->layerbits());
  return layeroutput[2];
}

unsigned long int HGCalCluster::CoreShowerLen()  {
  std::vector<int> layeroutput = ShowerLengthProperties(this->layerbits());
  return layeroutput[3];
}

unsigned long int HGCalCluster::E_EM_over_E_Fraction()  { 
  double E_EM_over_E_temp = Energy_ratio( this->e_em() , this->e() );
  double intpart;
  return modf(E_EM_over_E_temp,&intpart)*pow(2,8);
}

unsigned long int HGCalCluster::E_EM_over_E_Quotient()  { 
  double E_EM_over_E_temp = Energy_ratio( this->e_em() , this->e() );
  double intpart;
  double frac = modf(E_EM_over_E_temp,&intpart)*pow(2,8);
  return intpart;
}

unsigned long int HGCalCluster::E_EM_core_over_E_EM_Fraction()  { 
  double E_EM_core_over_E_EM_temp = Energy_ratio( this->e_em_core() , this->e() );
  double intpart;
  return modf(E_EM_core_over_E_EM_temp,&intpart)*pow(2,8);
}

unsigned long int HGCalCluster::E_EM_core_over_E_EM_Quotient()  { 
  double E_EM_core_over_E_EM_temp = Energy_ratio( this->e_em_core() , this->e() );
  double intpart;
  double frac =  modf(E_EM_core_over_E_EM_temp,&intpart)*pow(2,8);
  return intpart;
}

unsigned long int HGCalCluster::E_H_early_over_E_Fraction()  { 
  double E_H_early_over_E_temp = Energy_ratio( this->e_h_early() , this->e() );
  double intpart;
  return  modf(E_H_early_over_E_temp,&intpart)*pow(2,8);
}

unsigned long int HGCalCluster::E_H_early_over_E_Quotient()  { 
  double E_H_early_over_E_temp = Energy_ratio( this->e_h_early() , this->e() );
  double intpart;
  double frac =  modf(E_H_early_over_E_temp,&intpart)*pow(2,8);
  return intpart;
}


const HGCalCluster& HGCalCluster::operator+=(const HGCalCluster& c) {
  // Not handling field widths
  HGCalCluster original( *this );
  this->set_n_tc( this->n_tc() + c.n_tc() );
  this->set_e( this->e() + c.e() );
  this->set_e_em( this->e_em() + c.e_em() );
  this->set_e_em_core( this->e_em_core() + c.e_em_core() );
  this->set_e_h_early( this->e_h_early() + c.e_h_early() );
  this->set_w( this->w() + c.w() );
  this->set_n_tc_w( this->n_tc_w() + c.n_tc_w() );
  this->set_w2( this->w2() + c.w2() );
  this->set_wz( this->wz() + c.wz() );
  this->set_weta( this->weta() + c.weta() );
  this->set_wphi( this->wphi() + c.wphi() );
  this->set_wroz( this->wroz() + c.wroz() );
  this->set_wz2( this->wz2() + c.wz2() );
  this->set_weta2( this->weta2() + c.weta2() );
  this->set_wphi2( this->wphi2() + c.wphi2() );
  this->set_wroz2( this->wroz2() + c.wroz2() );

  this->set_layerbits( this->layerbits() | c.layerbits() );
  this->set_sat_tc( this->sat_tc() | c.sat_tc() );
  this->set_shapeq( this->shapeq() | c.shapeq() );

  if ( w_ <= 52438 && original.shapeq() == 1 && c.shapeq() == 1 ) { // Magic numbers
    this->set_shapeq(1);
  }
  else {
    this->set_shapeq(0);

    if ( this->w() > c.w() ) {
      this->set_w( original.w() );
      this->set_w2( original.w2() );
      this->set_wz( original.wz() );
      this->set_weta( original.weta() );
      this->set_wphi( original.wphi() );
      this->set_wroz( original.wroz() );
      this->set_wz2( original.wz2() );
      this->set_weta2( original.weta2() );
      this->set_wphi2( original.wphi2() );
      this->set_wroz2( original.wroz2() );
      this->set_n_tc_w( original.n_tc_w() );
    }
    else {
      this->set_w( c.w() );
      this->set_w2( c.w2() );
      this->set_wz( c.wz() );
      this->set_weta( c.weta() );
      this->set_wphi( c.wphi() );
      this->set_wroz( c.wroz() );
      this->set_wz2( c.wz2() );
      this->set_weta2( c.weta2() );
      this->set_wphi2( c.wphi2() );
      this->set_wroz2( c.wroz2() );
      this->set_n_tc_w( c.n_tc_w() );
    }

  }

  return *this;
}
