#ifndef L1Trigger_Phase2L1ParticleFlow_HTMHT_h
#define L1Trigger_Phase2L1ParticleFlow_HTMHT_h

#include "DataFormats/L1TParticleFlow/interface/jets.h"
#include "DataFormats/L1TParticleFlow/interface/sums.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/dbgPrintf.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/jetmet/L1SeedConePFJetEmulator.h"

#ifndef CMSSW_GIT_HASH
#include "hls_math.h"
#endif

#include <vector>
#include <numeric>
#include <algorithm>
#include "ap_int.h"
#include "ap_fixed.h"

namespace P2L1HTMHTEmu {
  typedef l1ct::pt_t pt_t;          // Type for pt/ht 1 unit = 0.25 GeV; max = 16 TeV
  typedef l1ct::glbeta_t etaphi_t;  // Type for eta & phi

  typedef ap_fixed<12, 3> radians_t;
  typedef ap_fixed<9, 2> cossin_t;
  typedef ap_fixed<16, 13> pxy_t;
  static constexpr int Fin = 6;  // Number of decimal bits/precision of met squared i.e. 2*16 - 2*13
  static constexpr int Fout = pt_t::width - pt_t::iwidth;  // Number of decimal bits/precision of output met

  static constexpr int N_TABLE = 2048;

  // Class for intermediate variables
  class PtPxPy {
  public:
    pt_t pt = 0.;
    pxy_t px = 0.;
    pxy_t py = 0.;

    PtPxPy operator+(const PtPxPy& b) const {
      PtPxPy c;
      c.pt = this->pt + b.pt;
      c.px = this->px + b.px;
      c.py = this->py + b.py;
      return c;
    }
  };

  namespace Scales {
    const ap_fixed<12, -4> scale_degToRad = M_PI / 180.;
  };  // namespace Scales

  template <class data_T, class table_T, int N>
  void init_sinphi_table(table_T table_out[N]) {
    for (int i = 0; i < N; i++) {
      double x = i * (M_PI / 180.) / 2.;
      table_T sin_x = std::sin(x);
      table_out[i] = sin_x;
    }
  }
  template <class in_t, class table_t, int N>
  table_t sine_with_conversion(etaphi_t hwPhi) {
    table_t sin_table[N];
    init_sinphi_table<in_t, table_t, N>(sin_table);
    table_t out = sin_table[hwPhi];
    return out;
  }

  // Software emulation of hls::atan2(pxy_t, pxy_t) = generic_atan2<W=16,I=13>.
  // Replicates the fixed-point CORDIC arithmetic bit-exactly so that the CMSSW
  // emulator matches the HLS firmware/csim output.
  inline ap_fixed<12, 3> atan2_cordic(pxy_t in1, pxy_t in2) {
    // Widths match generic_atan2<W=16,I=13>
    static constexpr int W  = pxy_t::width;   // 16
    static constexpr int I  = pxy_t::iwidth;  // 13
    static constexpr int WC = W + 7;          // 23  (CORDIC working width)
    static constexpr int NITER = WC - 3;      // 20  (fractional bits of CORDIC accumulator)

    // Constants — same types as generic_atan2
    static const ap_fixed<W + 1, 3> pi_ap("0x3.243F6A8885A308D3");
    static const ap_fixed<W + 7, 3> pi2_ap("0x1.921FB54442D1846");  // widened to WC,3
    static const ap_fixed<W + 1, 3> pi4_ap("0x0.C90FDAA22168C23");
    static const ap_fixed<W + 1, 3> pi3n_ap("-0x2.5B2F8FE6643A469");

    // CORDIC atan LUT in ap_fixed<WC,3> — populated once (AP_TRN default)
    static ap_fixed<WC, 3> atan_lut[NITER];
    static bool lut_ready = false;
    if (!lut_ready) {
      for (int i = 0; i < NITER; i++)
        atan_lut[i] = std::atan(std::pow(2.0, (double)-i));
      lut_ready = true;
    }

    ap_uint<2> signin1 = (in1 > 0) ? 2u : (in1 == 0) ? 1u : 0u;
    ap_uint<2> signin2 = (in2 > 0) ? 2u : (in2 == 0) ? 1u : 0u;

    ap_fixed<W, 3> out;

    // Corner cases (mirror generic_atan2 exactly)
    if (signin1 == 1 && signin2 == 2) { out = 0;        return out; }
    if (signin1 == 1 && signin2 == 0) { out = pi_ap;    return out; }
    if (signin1 == 2 && signin2 == 1) { out = pi2_ap;   return out; }
    if (signin1 == 0 && signin2 == 1) { out = -pi2_ap;  return out; }
    if (in1 == in2) {
      if      (signin1 == 2) { out = pi4_ap;  return out; }
      else if (signin1 == 1) { out = 0;       return out; }
      else                   { out = pi3n_ap; return out; }
    }

    // Absolute values: ap_fixed<W+1, I+1>
    ap_fixed<W + 1, I + 1> in1abs = (signin1 == 0) ? (ap_fixed<W + 1, I + 1>)(-in1) : (ap_fixed<W + 1, I + 1>)(in1);
    ap_fixed<W + 1, I + 1> in2abs = (signin2 == 0) ? (ap_fixed<W + 1, I + 1>)(-in2) : (ap_fixed<W + 1, I + 1>)(in2);

    // Bit-reinterpretation: ap_fixed<W+1,I+1> → ap_fixed<W+1,2> (mirrors the sft loop in generic_atan2)
    ap_fixed<W + 1, 2> in1abs_sft, in2abs_sft;
    in1abs_sft.range() = in1abs.range();
    in2abs_sft.range() = in2abs.range();

    ap_fixed<WC, 3> cx, cy, cz;
    if (in1abs > in2abs) { cx = in1abs_sft; cy = in2abs_sft; }
    else                 { cx = in2abs_sft; cy = in1abs_sft; }
    cz = 0;

    for (int i = 0; i < NITER; i++) {
      ap_fixed<WC, 3> cx_new, cy_new, cz_new;
      if (cy >= 0) {
        cx_new = cx + (cy >> i); cy_new = cy - (cx >> i); cz_new = cz + atan_lut[i];
      } else {
        cx_new = cx - (cy >> i); cy_new = cy + (cx >> i); cz_new = cz - atan_lut[i];
      }
      cx = cx_new; cy = cy_new; cz = cz_new;
    }

    if (in1abs > in2abs)
      cz = (ap_fixed<WC, 3>)pi2_ap - cz;

    if      (signin2 == 0 && signin1 == 2) out = (ap_fixed<WC, 3>)pi_ap - cz;
    else if (signin2 == 0 && signin1 == 0) out = cz - (ap_fixed<WC, 3>)pi_ap;
    else if (signin2 == 2 && signin1 == 0) out = -cz;
    else                                    out = cz;

    return out;
  }

  inline etaphi_t phi_cordic(pxy_t y, pxy_t x) {
#ifdef CMSSW_GIT_HASH
    ap_fixed<12, 3> phi = atan2_cordic(y, x);
#else
    ap_fixed<12, 3> phi = hls::atan2(y, x);
#endif
    ap_fixed<16, 9> etaphiscale = (float)l1ct::Scales::INTPHI_PI / M_PI;  // radians to hwPhi
    return phi * etaphiscale;
  }

  inline PtPxPy mht_compute(l1ct::Jet jet) {
    // Add an extra bit to px/py for the sign, and one additional bit to improve precision (pt_t is ap_ufixed<14, 12>)
    PtPxPy v_pxpy;

    //Initialize table once
    cossin_t sin_table[N_TABLE];
    init_sinphi_table<etaphi_t, cossin_t, N_TABLE>(sin_table);

    cossin_t sinphi;
    cossin_t cosphi;
    bool sign = jet.hwPhi.sign();

    etaphi_t hwphi = jet.hwPhi;

    // Reduce precision of hwPhi
    ap_int<10> phi;
    phi.V = hwphi(11, 1);
    phi = (phi > 0) ? phi : (ap_int<10>)-phi;  //Only store values for positive phi, pick up sign later

    sinphi = sin_table[phi];

    sinphi = (sign > 0) ? (cossin_t)(-sign * sinphi) : sinphi;  // Change sign bit if hwPt is negative, sin(-x)=-sin(x)
    cosphi = sin_table[phi + 90 * 2];  //cos(x)=sin(x+90). Do nothing with sign, cos(-θ) = cos θ,

    v_pxpy.pt = jet.hwPt;
    v_pxpy.py = jet.hwPt * sinphi;
    v_pxpy.px = jet.hwPt * cosphi;

    return v_pxpy;
  }
}  // namespace P2L1HTMHTEmu

//TODO replace with l1ct::Jet
inline l1ct::Sum htmht(std::vector<l1ct::Jet> jets) {
  // compute jet px, py
  std::vector<P2L1HTMHTEmu::PtPxPy> ptpxpy;
  ptpxpy.resize(jets.size());
  std::transform(
      jets.begin(), jets.end(), ptpxpy.begin(), [](const l1ct::Jet& jet) { return P2L1HTMHTEmu::mht_compute(jet); });

  // Sum pt, px, py over jets
  P2L1HTMHTEmu::PtPxPy hthxhy = std::accumulate(ptpxpy.begin(), ptpxpy.end(), P2L1HTMHTEmu::PtPxPy());

  // Compute the MHT magnitude and direction
  l1ct::Sum ht;
  ht.hwSumPt = hthxhy.pt;
#ifdef CMSSW_GIT_HASH
  double d = std::sqrt(((hthxhy.px * hthxhy.px) + (hthxhy.py * hthxhy.py)).to_double());
  // emulate hls::sqrt internal rounding
  double rounded = std::round(d * (1 << P2L1HTMHTEmu::Fin)) / (1 << P2L1HTMHTEmu::Fin);
  // emulate AP_TRN conversion to output type
  double truncated = std::floor(rounded * (1 << P2L1HTMHTEmu::Fout)) / (1 << P2L1HTMHTEmu::Fout);
  P2L1HTMHTEmu::pt_t hwPt_hls = truncated;
  ht.hwPt = hwPt_hls;
#else
  ht.hwPt = hls::sqrt(((hthxhy.px * hthxhy.px) + (hthxhy.py * hthxhy.py)));
#endif
  ht.hwPhi = P2L1HTMHTEmu::phi_cordic(hthxhy.py, hthxhy.px);
  return ht;
}

#endif
