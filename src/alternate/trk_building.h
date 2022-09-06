#ifndef __PROMPT_TRKBUILD_TRK_BUILDING_H__
#define __PROMPT_TRKBUILD_TRK_BUILDING_H__

//#include "../emtf_hlslib/layer_helpers.h"
#include <iostream>
#include "common.h"
#include "types.h"
#include "patterns.h"

namespace prompt_trkbuild {


// Table showing the first 36 features sent to NN
//
// feat       | ME1/1 | ME1/2 |  ME2  |  ME3  |  ME4  |  RE1  |  RE2  |  RE3  |  RE4  | GE1/1 | GE2/1 |  ME0
// -----------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------
// emtf_phi   |   *   |   *   |   *   |   *   |   *   |   *   |   *   |   *   |   *   |   *   |   *   |   *
// emtf_theta |   *   |   *   |   *   |   *   |   *   |   *   |   *   |   *   |   *   |   *   |   *   |   *
// emtf_bend  |   *   |   *   |   *   |   *   |   *   |       |       |       |       |       |       |   *
// emtf_qual  |   *   |   *   |   *   |   *   |   *   |       |       |       |       |       |       |   *
// emtf_time  |       |       |       |       |       |       |       |       |       |       |       |
//
// 4 additional features are: ph_median, th_median, trk_qual, trk_bx


template <typename T>
void median_of_three(const T a0, const T a1, const T a2, T& out) {
  // Takes median, or min if any value is invalid (255)

  // Stage 1: compare-swap if (wire[i] > wire[j]) swap(wire[j], wire[i])
  const T tmp_1_0 = (a0 > a1) ? a1 : a0;    // min(a0,a1) - def: a0
  const T tmp_1_1 = (a0 > a1) ? a0 : a1;    // max(a0,a1) - def: a1

  // Stage 2
  const T tmp_2_1 = (tmp_1_1 > a2) ? a2 : tmp_1_1; // min(a2, max(a0,a1))

  // Stage 3
  const T tmp_3_0 = (tmp_1_0 > tmp_2_1) ? tmp_2_1 : tmp_1_0;    // min
  const T tmp_3_1 = (tmp_1_0 > tmp_2_1) ? tmp_1_0 : tmp_2_1;    // middle

  // Output - if any are invalid, take the min
  const T invalid_val = -1;
  out = (a0==invalid_val or a1==invalid_val or a2==invalid_val) ? tmp_3_0 : tmp_3_1;
}



void trk_building( best_trk_t best_trks[num_tracks], const emtf_seg_t emtf_segs[num_sites][max_ch_site][max_seg_ch],
                   emtf_seg_t trk_features[num_tracks][max_num_layers], 
                   emtf_phi_t ph_median_out[num_tracks], emtf_theta_t th_median_out[num_tracks]){
    
    for(int trk=0; trk<num_tracks; trk++){
    
        trk_zone_t trk_zone = best_trks[trk].zone;
        trk_qual_t trk_qual = best_trks[trk].qual;
        trk_patt_t trk_patt = best_trks[trk].patt;
        trk_col_t  trk_col = best_trks[trk].col;

        const int hitmap_offset = 27;
        const int bits_to_shift = 4;
        const int middle_of_range = (1 << (bits_to_shift-1)); // 8

        // Get Phi Median from pattern col
        trk_col_t col_corr = trk_col + hitmap_offset;  // 27 bc of start of col offset
        emtf_phi_t ph_median = (static_cast<emtf_phi_t>(col_corr) << bits_to_shift) + middle_of_range;

        ph_median_out[trk] = ph_median;
        

        if(trk_qual > 0){

            emtf_seg_t best_phi_hits[num_sites];
            emtf_phi_t best_phi_hits_diffs[num_sites];

            for(int site=0; site<num_sites; site++)
                best_phi_hits[site].seg_valid = 0;
        

            /****** Add hits using Phi only *********/

            // Add Hits Based on proximity to center of pattern in Phi (and verify that it is within pattern)
            for(int row=0; row<num_hitmap_rows; row++){

                // Get pattern columns with padding offset and trk_col offset
                trk_col_t  patt_start_col_padded = col_corr + col_start_table[trk_zone][trk_patt][row];
                trk_col_t  patt_mid_col_padded = col_corr + col_mid_table[trk_zone][trk_patt][row];
                trk_col_t  patt_stop_col_padded = col_corr + col_stop_table[trk_zone][trk_patt][row];

                // Get pattern middle with trk_col offset, but without padding.
                trk_col_t patt_mid_col = (patt_mid_col_padded <= max_padding)? ap_uint<9>(0) : ap_uint<9>(patt_mid_col_padded - max_padding);

                // Calculate middle of pattern range (for this zone-patt-row) at phi resolution
                emtf_phi_t patt_mid_ph = (static_cast<emtf_phi_t>(patt_mid_col) << bits_to_shift) + middle_of_range;

                // Loop through all sites in a row (they match separately)
                for(int site=0; site<max_sites_in_row; site++){

                    int curr_site = sites_by_row[trk_zone][row][site];
                    emtf_phi_t min_diff = -1;
                    ap_uint<1> min_diff_sign = 0;
                    emtf_seg_t min_seg;
                    min_seg.seg_valid = 0; // incase no valid segments

                    if(curr_site != no_site){
                        // Loop through the segments and select the best one (this is actually done in parallel)
                        for(int ch=0; ch<max_ch_site; ch++){
                            for(int iseg=0; iseg<max_seg_ch; iseg++){

                                emtf_seg_t curr_seg = emtf_segs[curr_site][ch][iseg];
                                emtf_phi_t diff;
                                ap_uint<1> diff_sign;

                                if(curr_seg.seg_valid and curr_seg.seg_zones[num_zones-1-trk_zone] == 1){
                                    // abs diff
                                    if(patt_mid_ph > curr_seg.emtf_phi){
                                        diff = patt_mid_ph - curr_seg.emtf_phi;
                                    }
                                    else{
                                        diff = curr_seg.emtf_phi - patt_mid_ph;
                                    }

                                    if(diff < min_diff){
                                        min_diff = diff;
                                        min_seg = curr_seg;
                                    }
                                }
                            }
                        }

                        // Check that segment is inside of pattern
                        trk_col_t min_seg_col_padded = (min_seg.emtf_phi >> 4) + max_padding; //real 315bit col w/ padding

                        if(min_seg.seg_valid and min_seg_col_padded >= patt_start_col_padded and min_seg_col_padded <= patt_stop_col_padded){
                            // match is valid
                            best_phi_hits[curr_site] = min_seg;
                            best_phi_hits_diffs[curr_site] = min_diff;
                        }
                    }
                } // site
            } // row

            
            // Layer order
            // | ME1/1 | ME1/2 |  ME2  |  ME3  |  ME4  |  RE1  |  RE2  |  RE3  |  RE4  | GE1/1 | GE2/1 |  ME0

            // Reduce sites to layers and reorder them to Jia Fu's layer order
            for(int layer=0; layer<max_num_layers; layer++){
                if(map_sites_to_layers[layer][1] == no_site){
                    int site = map_sites_to_layers[layer][0];
                    trk_features[trk][layer] = best_phi_hits[site];
                }
                else{
                    // Handle merge of 2 sites into 1 layer
                    int site0 = map_sites_to_layers[layer][0];
                    int site1 = map_sites_to_layers[layer][1];

                    if( best_phi_hits[site0].seg_valid and !best_phi_hits[site1].seg_valid){
                        trk_features[trk][layer] = best_phi_hits[site0]; // only site0 is valid
                    }else if( best_phi_hits[site1].seg_valid and !best_phi_hits[site0].seg_valid){
                        trk_features[trk][layer] = best_phi_hits[site1]; // only site1 is valid
                    }else if( best_phi_hits[site1].seg_valid and best_phi_hits[site0].seg_valid){
                        // both sites are valid, pick smallest diff
                        if( best_phi_hits_diffs[site0] < best_phi_hits_diffs[site1])
                            trk_features[trk][layer] = best_phi_hits[site0]; // Site 0 is better or equal
                        else
                            trk_features[trk][layer] = best_phi_hits[site1]; // Site 1 is better
                    }else{
                        trk_features[trk][layer].seg_valid = 0; // neither site is valid, no attached hit
                    }
                }
            }
            

            /************************************* Find Theta 'Median' **********************************************/

            constexpr const int trk_theta_indices[3][9] = {
                                    { 2,  3,  4, 2, 3, 4, 10,  7,  8}, // [ME2_t1, ME3_t1, ME4_t1, ME2_t2, ME3_t2, ME4_t2, GE2, RE3, RE4]
                                    { 2,  3,  4, 2, 3, 4,  6,  7,  8}, // [ME2_t1, ME3_t1, ME4_t1, ME2_t2, ME3_t2, ME4_t2, RE2, RE3, RE4]
                                    { 1,  0, 11, 1, 0, 11,  5,  9, 11} // [ME12_t1, ME11_t1, ME0_t2, ME12_t2, ME11_t2, ME0_t2, RE1, GE1, ME0_t1]
            };

            constexpr const int trk_theta_which_th[3][9]     = { {1,1,1, 2,2,2, 1,1,1}, {1,1,1, 2,2,2, 1,1,1}, {1,1,2, 2,2,2, 1,1,1} };

            const emtf_theta_t invalid_median = -1; // invalid for sorting, set all bits to 1

            emtf_theta_t thetas[3][9];
            for(int i=0; i<3; i++){
                for(int j=0; j<9; j++){
                    int layer = trk_theta_indices[i][j];
                    int which_th = trk_theta_which_th[i][j];

                    if(trk_features[trk][layer].seg_valid)
                        thetas[i][j] = (which_th == 1)? trk_features[trk][layer].emtf_theta1 : trk_features[trk][layer].emtf_theta2; 
                    
                    // all invalid numbers get set to highest number possible (255 for 8 bit)
                    if(!trk_features[trk][layer].seg_valid or thetas[i][j] == 0)
                        thetas[i][j] = invalid_median; // largest number possible
                }
            }

            // For each option (non-st1 default, non-st1 alternative, st1), Find median of every 3, then find median of those 
            emtf_theta_t median_options[3];
            for(int i=0; i<3; i++){
                emtf_theta_t stage1_1, stage1_2, stage1_3, stage2_1;
                median_of_three(thetas[i][0],thetas[i][1],thetas[i][2], stage1_1);
                median_of_three(thetas[i][3],thetas[i][4],thetas[i][5], stage1_2);
                median_of_three(thetas[i][6],thetas[i][7],thetas[i][8], stage1_3);
                median_of_three(stage1_1, stage1_2, stage1_3, median_options[i]);
            }
            
            emtf_theta_t th_median;
            if( trk_zone != 2)
                th_median = median_options[0];
            else
                th_median = median_options[1];

            // If above fails, try station 1 median
            if(th_median == invalid_median)
                th_median = median_options[2];
            
            // If all of theta fails, the median is invalid
            if(th_median == invalid_median)
                th_median = 0; // invalid theta value
                
            th_median_out[trk] = th_median;
            


            /************************************ Compare Thetas to  th_median *************************/

            const int th_window = 8; // if diff <= th_window, it is invalid

            for(int layer=0; layer<max_num_layers; layer++){
                emtf_seg_t curr_hit = trk_features[trk][layer];

                if(curr_hit.seg_valid){
                    // set with invalid values
                    emtf_theta_t diff1 = th_window + 1;
                    emtf_theta_t diff2 = th_window + 1;

                    // Calc abs theta 1 diff
                    if(curr_hit.emtf_theta1 != 0){
                        if(th_median > curr_hit.emtf_theta1)
                            diff1 = th_median - curr_hit.emtf_theta1;
                        else
                            diff1 = curr_hit.emtf_theta1 - th_median;
                    }

                    // Calc abs theta 2 diff
                    if(curr_hit.emtf_theta2 != 0){
                        if(th_median > curr_hit.emtf_theta2)
                            diff2 = th_median - curr_hit.emtf_theta2;
                        else
                            diff2 = curr_hit.emtf_theta2 - th_median;
                    }

                    if(diff2 < diff1 and diff2 < th_window){
                        // Attach theta 2 as correct theta
                        trk_features[trk][layer].emtf_theta1 = trk_features[trk][layer].emtf_theta2;
                    }
                    else if( diff1 >= th_window){
                        trk_features[trk][layer].emtf_theta1 = 0; // invalid
                        trk_features[trk][layer].seg_valid = 0;
                    }
                }
            }

        } // if qual valid
        else{
            th_median_out[trk] = 0;
            ph_median_out[trk] = 0;
            for(int layer=0; layer<max_num_layers; layer++)
                trk_features[trk][layer].seg_valid = 0;
        }
    } // trk


}


} // namespace


#endif
