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
                   trk_feat_t trk_features[num_tracks][num_emtf_features], seg_id_t trk_seg_ids[num_tracks][max_num_layers]){
                   
    // Zero the outputs
    for(int trk=0; trk<num_tracks; trk++){
        for(int ifeat=0; ifeat<num_emtf_features; ifeat++)
            trk_features[trk][ifeat] = 0;

        for(int layer=0; layer<max_num_layers; layer++){
            trk_seg_ids[trk][layer].isite = 0;
            trk_seg_ids[trk][layer].ichamber = 0;
            trk_seg_ids[trk][layer].iseg = 0;
            trk_seg_ids[trk][layer].ibx = 0;
        }
    }
        
    
    // Process each track
    for(int trk=0; trk<num_tracks; trk++){
    
        trk_zone_t trk_zone = best_trks[trk].zone;
        trk_qual_t trk_qual = best_trks[trk].qual;
        trk_patt_t trk_patt = best_trks[trk].patt;
        trk_col_t  trk_col = best_trks[trk].col;

        emtf_phi_t ph_median = 0;

        // Dont need to run for invalid trk
        if(trk_qual == 0)
            continue;

        const int hitmap_offset = 27;
        const int bits_to_shift = 4;
        const int middle_of_range = (1 << (bits_to_shift-1)); // 8

        /****** Get Phi Median from trk col *******/
        trk_col_t col_corr = trk_col + hitmap_offset;  // 27 bc of start of col offset
        ph_median = (static_cast<emtf_phi_t>(col_corr) << bits_to_shift) + middle_of_range;

        /***** Get the pattern info for each row ****/
        trk_col_t  patt_start_col_padded[num_hitmap_rows]; 
        trk_col_t  patt_mid_col_padded[num_hitmap_rows]; 
        trk_col_t  patt_stop_col_padded[num_hitmap_rows]; 
        trk_col_t  patt_mid_col[num_hitmap_rows];
        emtf_phi_t patt_mid_ph[num_hitmap_rows];

        for(int row=0; row<num_hitmap_rows; row++){
            // Get pattern columns with padding offset and trk_col offset
            patt_start_col_padded[row] = col_corr + col_start_table[trk_zone][trk_patt][row];
            patt_mid_col_padded[row] = col_corr + col_mid_table[trk_zone][trk_patt][row];
            patt_stop_col_padded[row] = col_corr + col_stop_table[trk_zone][trk_patt][row];

            // Get column in middle of pattern rows for this tracks zone (for full 315col sector)
            patt_mid_col[row] = (patt_mid_col_padded[row] <= max_padding)? ap_uint<9>(0) : ap_uint<9>(patt_mid_col_padded[row] - max_padding);

            // Calculate middle of pattern range (for this zone-patt-row) at phi resolution
            patt_mid_ph[row] = (static_cast<emtf_phi_t>(patt_mid_col[row]) << bits_to_shift) + middle_of_range;
        }


        /******* Updated method with chamber gating ********/
        const int site_row_by_zone[num_zones][num_sites] = {
            // ME0 GE11 ME11 ME12 RE12 GE21 RE22 ME21 ME22 ME31 ME32 RE31 RE32 ME41 ME42 RE41 RE42
            {  0,  1,   2,   -1,  -1,  3,   -1,  4,   -1,  5,   -1,  6,   -1,  7,   -1,  7,   -1}, // Zone 0
            {  -1, 0,   1,   2,   2,   3,   -1,  4,   -1,  5,   5,   6,   6,   7,   7,   7,   7}, // Zone 1
            {  -1, -1,  -1,  0,   1,   -1,  2,   -1,  3,   -1,  4,   -1,  5,   -1,  6,   -1,  7}, // Zone 2
        };


        int max_ch_per_gate = 3;
        int ch_per_gate_bt[num_site_types] = {3,2,2}; // Ten degree chambers have max of 3 possible chambers
        const int num_gates = 8; 

        const int gate_edges[num_gates-1] = {101,120,158,176,195,233,251}; // where chamber gate changes

        const int possible_chambers_start_bt[num_site_types][num_gates] = {
            //edge: 0   1   2   3   4   5   6   
            //      |   |   |   |   |   |   |
                {0,   0,  1,  2,  2,  3,  4,  4}, // ten degree chamber start
                {0,   0,  1,  1,  1,  2,  2,  2}, // twenty degree chamber start
                {0,   1,  1,  1,  2,  2,  2,  3}  // ME0 chambers
        };
        



        /********* Attach Best Segment for each layer *********/


        emtf_seg_t min_seg[num_sites];
        emtf_phi_t min_diff[num_sites]; 
        int min_seg_site_idx[num_sites];
        int min_seg_ch_idx[num_sites];
        int min_seg_seg_idx[num_sites];

        // Get the best segment for each site
        for(int site=0; site<num_sites; site++){
            min_seg[site].seg_valid = 0; 
            min_diff[site] = static_cast<emtf_phi_t>(-1); // set diff to max

            // Get the row for this site
            int site_row = site_row_by_zone[trk_zone][site];
            if(site_row == -1) // unused site in this zone
                continue;

            // Get the current gate and the possible chambers with valid segments
            int curr_gate = 0;
            for(int i=0; i<num_gates-1; i++){
                if(patt_mid_col[site_row] >= gate_edges[i])
                    curr_gate = i+1; 
            }

            int site_type = type_by_site[site];
            int start_ch = possible_chambers_start_bt[site_type][curr_gate]; // Gets starting chamber to look

            // Find the best match for this site
            for(int ch_in_gate=0; ch_in_gate<ch_per_gate_bt[site_type]; ch_in_gate++){

                // The real chamber id in current site
                int curr_ch = start_ch + ch_in_gate;

                for(int iseg=0; iseg<max_seg_ch; iseg++){

                    emtf_seg_t curr_seg = emtf_segs[site][curr_ch][iseg];

                    // Check that segment is inside of pattern
                    trk_col_t curr_seg_col_padded = (curr_seg.emtf_phi >> 4) + max_padding; //real 315bit col w/ padding

                    bool good_segment = (curr_seg.seg_valid and                                    // Segment is valid 
                                        curr_seg.seg_zones[num_zones-1-trk_zone] == 1 and          // Segment is within this zone
                                        curr_seg_col_padded >= patt_start_col_padded[site_row] and // Segment is within pattern
                                        curr_seg_col_padded <= patt_stop_col_padded[site_row]);

                    if(good_segment){
                        // abs diff
                        emtf_phi_t diff;
                        if(patt_mid_ph[site_row] > curr_seg.emtf_phi)
                            diff = patt_mid_ph[site_row]- curr_seg.emtf_phi;
                        else
                            diff = curr_seg.emtf_phi - patt_mid_ph[site_row];

                        if(diff < min_diff[site]){
                            min_diff[site] = diff;
                            min_seg[site] = curr_seg;
                            min_seg_site_idx[site] = site;
                            min_seg_ch_idx[site] = curr_ch;
                            min_seg_seg_idx[site] = iseg;
                        }
                    }
                } // iseg
            } // chamber
        } // Site


        // Some sites are combined into 1 layer for the NN (ordered based on Jia Fus NN)
        const bool multisite_layers[max_num_layers] = {false, false, true, true, true, false, false, true, true, false, false, false};

        // Hold the attached segments
        attached_seg_t attached_segs[max_num_layers];

        for(int layer=0; layer<max_num_layers; layer++){

            int site0 = map_layers_to_sites[layer][0];
            int site1 = map_layers_to_sites[layer][1]; // may be none (-1)
            int best_site;

            if(site1 == no_site)
                best_site = site0;
            else{
                int priority_site = ((min_seg_ch_idx[site0] * 2) <= (min_seg_ch_idx[site1]+1) )? site0 : site1; // If equal, this one gets priority
                int low_prio_site = ((min_seg_ch_idx[site0] * 2) <= (min_seg_ch_idx[site1]+1) )? site1 : site0;
                 
                if(min_diff[priority_site] <= min_diff[low_prio_site])
                    best_site = priority_site;
                else
                    best_site = low_prio_site;
            }

            attached_segs[layer].emtf_seg = min_seg[best_site];
            attached_segs[layer].seg_id.isite = min_seg_site_idx[best_site];
            attached_segs[layer].seg_id.ichamber = min_seg_ch_idx[best_site];
            attached_segs[layer].seg_id.iseg = min_seg_seg_idx[best_site];
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

                if(attached_segs[layer].emtf_seg.seg_valid)
                    thetas[i][j] = (which_th == 1)? attached_segs[layer].emtf_seg.emtf_theta1 : attached_segs[layer].emtf_seg.emtf_theta2; 
                
                // all invalid numbers get set to highest number possible (255 for 8 bit)
                if(!attached_segs[layer].emtf_seg.seg_valid or thetas[i][j] == 0)
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
            

        /************************************ Compare Thetas to  th_median *************************/
        const int th_window = 8; // if diff <= th_window, it is invalid

        for(int layer=0; layer<max_num_layers; layer++){
            emtf_seg_t curr_hit = attached_segs[layer].emtf_seg;

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
                    attached_segs[layer].emtf_seg.emtf_theta1 = attached_segs[layer].emtf_seg.emtf_theta2;
                }
                else if( diff1 >= th_window){
                    attached_segs[layer].emtf_seg.emtf_theta1 = 0; // invalid
                    attached_segs[layer].emtf_seg.seg_valid = 0;
                }
            }
        }


        /*******  Prepare the Track Features Output ***********/

        // Collect Segment Features Phi, Theta, and valid bits
        for(int layer=0; layer<max_num_layers; layer++){
            if(attached_segs[layer].emtf_seg.seg_valid){
                trk_features[trk][layer] = static_cast<trk_feat_t>(attached_segs[layer].emtf_seg.emtf_phi) - ph_median; // get signed phi difference
                trk_features[trk][layer+max_num_layers] =  static_cast<trk_feat_t>(attached_segs[layer].emtf_seg.emtf_theta1) - th_median; // get signed theta difference
            }
        }
          
        // Collect the rest of the features manually for now
        trk_features[trk][24] = (attached_segs[0].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[0].emtf_seg.emtf_bend) : static_cast<trk_feat_t>(0); // Bends are sign extended
        trk_features[trk][25] = (attached_segs[1].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[1].emtf_seg.emtf_bend) : static_cast<trk_feat_t>(0);
        trk_features[trk][26] = (attached_segs[2].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[2].emtf_seg.emtf_bend) : static_cast<trk_feat_t>(0);
        trk_features[trk][27] = (attached_segs[3].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[3].emtf_seg.emtf_bend) : static_cast<trk_feat_t>(0);
        trk_features[trk][28] = (attached_segs[4].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[4].emtf_seg.emtf_bend) : static_cast<trk_feat_t>(0);
        trk_features[trk][29] = (attached_segs[11].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[11].emtf_seg.emtf_bend) : static_cast<trk_feat_t>(0);
        
        trk_features[trk][30] = (attached_segs[0].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[0].emtf_seg.emtf_qual) : static_cast<trk_feat_t>(0);
        trk_features[trk][31] = (attached_segs[1].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[1].emtf_seg.emtf_qual) : static_cast<trk_feat_t>(0);
        trk_features[trk][32] = (attached_segs[2].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[2].emtf_seg.emtf_qual) : static_cast<trk_feat_t>(0);
        trk_features[trk][33] = (attached_segs[3].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[3].emtf_seg.emtf_qual) : static_cast<trk_feat_t>(0);
        trk_features[trk][34] = (attached_segs[4].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[4].emtf_seg.emtf_qual) : static_cast<trk_feat_t>(0);
        trk_features[trk][35] = (attached_segs[11].emtf_seg.seg_valid)? static_cast<trk_feat_t>(attached_segs[11].emtf_seg.emtf_qual)  : static_cast<trk_feat_t>(0);
        
        trk_features[trk][36] = (static_cast<trk_feat_t>(ph_median) - 2744); // centers the ph_median range at 0 -> Note: 2744 = ((288cols/2 + 27) << 4) + 8
        trk_features[trk][37] = static_cast<trk_feat_t>(th_median);
        trk_features[trk][38] = static_cast<trk_feat_t>(trk_qual);
        trk_features[trk][39] = 0; //trk bx

    } // trk


}


} // namespace


#endif
