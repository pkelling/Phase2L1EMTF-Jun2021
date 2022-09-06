#ifndef __PROMPT_TRKBUILD_TYPES_H__
#define __PROMPT_TRKBUILD_TYPES_H__

// EMTF HLS
//#include "common.h"
//#include "traits.h"

namespace prompt_trkbuild {

// Segment datatypes
//
// Name          | Bit width | Is signed
// --------------|-----------|-----------
// emtf_phi_t    | 13        | 0
// emtf_bend_t   | 10        | 1
// emtf_theta1_t | 8         | 0
// emtf_theta2_t | 8         | 0
// emtf_qual1_t  | 4         | 0
// emtf_qual2_t  | 4         | 0
// emtf_time_t   | 4         | 1
// seg_zones_t   | 3         | 0
// seg_tzones_t  | 3         | 0
// seg_cscfr_t   | 1         | 0
// seg_gemdl_t   | 1         | 0
// seg_bx_t      | 2         | 1
// seg_valid_t   | 1         | 0
//
// Track datatypes
//
// Name          | Bit width | Is signed
// --------------|-----------|-----------
// trk_qual_t    | 6         | 0
// trk_patt_t    | 3         | 0
// trk_col_t     | 9         | 0
// trk_zone_t    | 2         | 0
// trk_tzone_t   | 2         | 0
// trk_gate_t    | 2         | 0
// trk_seg_t     | 8         | 0
// trk_seg_v_t   | 12        | 0
// trk_feat_t    | 13        | 1
// trk_invpt_t   | 14        | 1
// trk_phi_t     | 14        | 1
// trk_eta_t     | 14        | 1
// trk_d0_t      | 14        | 1
// trk_z0_t      | 14        | 1
// trk_beta_t    | 14        | 1
// trk_valid_t   | 1         | 0


// ******** Input Types (EMTF Segent) **********
typedef ap_uint<13> emtf_phi_t;
typedef ap_int<10>  emtf_bend_t;
typedef ap_uint<8>	emtf_theta_t;
typedef ap_uint<4>	emtf_qual_t;
typedef ap_int<4>	emtf_time_t;
typedef ap_uint<3>	seg_zones_t;
typedef ap_int<2>	seg_bx_t;
typedef ap_uint<1>	seg_valid_t;

typedef struct {
	emtf_phi_t 		emtf_phi;
	emtf_bend_t 	emtf_bend;
	emtf_theta_t 	emtf_theta1;
	emtf_theta_t 	emtf_theta2;
	emtf_qual_t 	emtf_qual;
	emtf_time_t 	emtf_time;
	seg_zones_t 	seg_zones;
	seg_bx_t 		seg_bx;
	seg_valid_t 	seg_valid;
} emtf_seg_t;


// ************ Output Data Types ******************
typedef ap_uint<6> trk_qual_t;
typedef ap_uint<3> trk_patt_t;
typedef ap_uint<9> trk_col_t;
typedef ap_uint<2> trk_zone_t;

/*
typedef ap_uint<2> trk_gate_t;
typedef ap_uint<8> trk_seg_t;
typedef ap_uint<12> trk_seg_v_t;
typedef ap_int<13> trk_feat_t;
typedef ap_int<14> trk_invpt_t;
typedef ap_int<14> trk_phi_t;
typedef ap_int<14> trk_eta_t;
typedef ap_int<14> trk_d0_t;
typedef ap_int<14> trk_z0_t;
typedef ap_int<14> trk_beta_t;
typedef ap_uint<1> trk_valid_t;

// Select the widest track datatype
typedef trk_invpt_t model_out_t;
*/


// Stage Types
typedef ap_uint<315> raw_hitmap_t;
typedef ap_uint<288> hitmap_t;

typedef struct {
    trk_patt_t patt;
    trk_qual_t qual;
} patt_match_t;

typedef struct {
    trk_patt_t patt;
    trk_qual_t qual;
    trk_col_t  col;
    trk_zone_t zone;
} best_trk_t;



} // namespace


#endif  // __EMTF_HLSLIB_TYPES_H__ not defined
